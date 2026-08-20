#!/usr/bin/env python
# coding: utf-8

# ## 1 · Imports and configuration

# In[ ]:


import io
import os
import re
import sys
import time
import json
import textwrap
from datetime import datetime, timezone

import pandas as pd
from Bio import Entrez, Medline


Entrez.email   = ""
Entrez.tool    = "ESCA-gene-prioritisation-litscreen"
Entrez.api_key = os.environ.get("NCBI_API_KEY")  


Entrez.max_tries          = 3
Entrez.sleep_between_tries = 15


RETMAX = 10        
SORT   = "relevance" 



UNSAFE_ALIAS_CHARS = set('()[]"\'{}&|:*~')

ALIAS_MIN_LEN = 4

ALIAS_BLOCKLIST = {
    "TS", "MS", "PC", "CAD", "APC", "NEU", "LOX", "CAT", "AGE", "ARM", "BAD",
    "CAP", "CAR", "COPE", "DDT", "FAT", "GAS", "HR", "ICE", "IMPACT", "LARGE",
    "MARS", "NELL", "PIGS", "REST", "SET", "SHE", "SPARC", "STAR", "TANK",
    "WARS", "CS", "PS", "AR", "ER", "PR", "MET", "SRC", "RAN", "ACHE", "MICE",
}

EXCLUDE_REVIEWS = True

INCLUDE_AMBIGUOUS_DISEASE_ACRONYMS = True

def build_disease_block(include_ambiguous=INCLUDE_AMBIGUOUS_DISEASE_ACRONYMS):
    """The esophageal-cancer half of the query. Unchanged in substance from
    your original, only restructured so it can be unit-tested independently."""
    terms = [
        '"Esophageal Neoplasms"[MeSH Terms]',
        '"esophageal cancer"[Title/Abstract]',
        '"oesophageal cancer"[Title/Abstract]',
        '"esophageal carcinoma"[Title/Abstract]',
        '"oesophageal carcinoma"[Title/Abstract]',
        '"esophageal adenocarcinoma"[Title/Abstract]',
        '"oesophageal adenocarcinoma"[Title/Abstract]',
        '"esophageal squamous cell carcinoma"[Title/Abstract]',
        '"oesophageal squamous cell carcinoma"[Title/Abstract]',
        'ESCC[Title/Abstract]',
    ]
    if include_ambiguous:
        terms += ['EAC[Title/Abstract]', 'ESCA[Title/Abstract]']
    return "(" + " OR ".join(terms) + ")"


# ## 2 · HGNC

# In[ ]:


def detect_multivalue_delimiter(series):
    joined = " ".join(series.dropna().astype(str).tolist()[:5000])
    n_pipe, n_comma = joined.count("|"), joined.count(",")
    return "|" if n_pipe > n_comma else ","


def split_multivalue(cell, delim):
    if not cell:
        return []
    return [tok.strip().strip('"').strip() for tok in cell.split(delim) if tok.strip()]


def load_hgnc(path):
    hgnc = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    required = {"Approved symbol", "Previous symbols", "Alias symbols"}
    missing = required - set(hgnc.columns)
    if missing:
        raise ValueError(f"custom.tsv is missing required column(s): {missing}")

    if "Status" in hgnc.columns:
        before = len(hgnc)
        hgnc = hgnc[hgnc["Status"].str.strip() == "Approved"].copy()
        print(f"  Status filter: kept {len(hgnc)}/{before} rows (Approved only)")
    else:
        print("  WARNING: no 'Status' column found — cannot filter withdrawn entries")

    delim = detect_multivalue_delimiter(
        pd.concat([hgnc["Previous symbols"], hgnc["Alias symbols"]])
    )
    print(f"  Multi-value delimiter detected: {delim!r}")
    return hgnc, delim


# ## 3 · Alias

# In[ ]:


def sanitize_alias(alias, approved_symbol, approved_symbols_set):
    a = alias.strip()

    if not a:
        return False, "empty"

    # Rule 0 — the gene's own approved symbol is exempt from all other rules.
    if a == approved_symbol:
        return True, "ok (approved symbol, exempt)"

    # Rule 1 — Boolean-significant characters. THE CRITICAL RULE.
    bad = sorted(set(a) & UNSAFE_ALIAS_CHARS)
    if bad:
        return False, f"unsafe char(s) {''.join(bad)} would break Boolean parsing"

    if a in approved_symbols_set:
        return False, "is the approved symbol of a different gene"

    # Rule 3 — too short to be specific.
    if len(a) < ALIAS_MIN_LEN:
        return False, f"shorter than ALIAS_MIN_LEN={ALIAS_MIN_LEN}"

    # Rule 4 — known ambiguous token / English word.
    if a.upper() in ALIAS_BLOCKLIST:
        return False, "on ALIAS_BLOCKLIST (ambiguous English word/abbreviation)"


    if "-" in a:
        return True, "ok (WARNING: hyphenated — may miss PubMed's phrase index)"

    return True, "ok"


def build_alias_dict(hgnc, delim):
    """Build {approved_symbol: [safe aliases]} plus a full audit trail."""
    approved_symbols = set(hgnc["Approved symbol"].str.strip())

    alias_dict, audit = {}, []

    for _, row in hgnc.iterrows():
        symbol = row["Approved symbol"].strip()

        candidates = [symbol]
        candidates += split_multivalue(row["Previous symbols"], delim)
        candidates += split_multivalue(row["Alias symbols"], delim)

        kept = []
        for a in dict.fromkeys(candidates):        # de-dupe, preserve order
            ok, reason = sanitize_alias(a, symbol, approved_symbols)
            audit.append({"gene": symbol, "alias": a,
                          "kept": ok, "reason": reason})
            if ok:
                kept.append(a)

        alias_dict[symbol] = sorted(set(kept))

    return alias_dict, pd.DataFrame(audit)


# ## 4 · Query construction and validation

# In[ ]:


def build_query(aliases, disease_block=None, exclude_reviews=EXCLUDE_REVIEWS):
    if disease_block is None:
        disease_block = build_disease_block()
    if not aliases:
        raise ValueError("alias list is empty — refusing to build a query")

    gene_block = "(" + " OR ".join(
        f'"{a}"[Title/Abstract]' for a in aliases
    ) + ")"

    q = f"{gene_block} AND {disease_block}"
    if exclude_reviews:
        q += " NOT review[Publication Type]" 
    return q


def build_query_ORIGINAL(aliases):
    gene_query = "(" + " OR ".join(
        f'"{alias}"[Title/Abstract]' for alias in aliases
    ) + ")"
    return (
        f'{gene_query} AND ('
        f'"Esophageal Neoplasms"[MeSH Terms] OR '
        f'"esophageal cancer"[Title/Abstract] OR '
        f'ESCC[Title/Abstract] OR ESCA[Title/Abstract]'
        f') AND NOT review[Publication Type]'
    )

QUOTED_PHRASE_RE = re.compile(r'"([^"]*)"')

def assert_query_wellformed(q):
    """Raise ValueError if the query is structurally unsound. Cheap and local."""
    problems = []

    depth = 0
    for ch in q:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                problems.append("closing ')' before a matching '('")
                break
    if depth != 0:
        problems.append(f"unbalanced parentheses (final depth {depth})")

    for phrase in QUOTED_PHRASE_RE.findall(q):
        bad = sorted(set(phrase) & UNSAFE_ALIAS_CHARS)
        if bad:
            problems.append(
                f'quoted phrase "{phrase}" contains grouping char(s) '
                f'{"".join(bad)} — will regroup the query if unquoted by PubMed'
            )

    if q.count('"') % 2 != 0:
        problems.append(f"odd number of double quotes ({q.count(chr(34))})")

    if re.search(r"\b(AND|OR|NOT)\s+(AND|OR|NOT)\b", q):
        problems.append("two Boolean operators adjacent (e.g. 'AND NOT')")

    if "()" in q.replace(" ", ""):
        problems.append("empty parenthesis group '()'")

    if problems:
        raise ValueError("Malformed query: " + "; ".join(problems) + f"\nQuery: {q}")
    return True


# ## 5 · Instrumented E-utilities

# In[ ]:


def esearch_checked(query, retmax=RETMAX, sort=SORT):
    """Run one ESearch and return a fully instrumented result dict."""
    assert_query_wellformed(query)                    

    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax, sort=sort)
    rec = Entrez.read(handle)
    handle.close()

    warn = rec.get("WarningList", {}) or {}
    err  = rec.get("ErrorList",   {}) or {}

    flags = []
    for key in ("PhraseNotFound", "FieldNotFound"):
        if err.get(key):
            flags.append(f"ERROR:{key}={list(err[key])}")
    for key in ("QuotedPhraseNotFound", "PhraseIgnored", "OutputMessage"):
        if warn.get(key):
            flags.append(f"WARN:{key}={list(warn[key])}")

    return {
        "count":             int(rec["Count"]),
        "pmids":             list(rec["IdList"]),
        "query_translation": rec.get("QueryTranslation", ""),
        "flags":             " | ".join(flags),
        "clean":             not flags,
    }


YEAR_RE = re.compile(r"\b(1[89]\d{2}|20\d{2})\b")

def extract_year(dp):
    m = YEAR_RE.search(dp or "")
    return m.group(1) if m else ""


def efetch_by_pmid(pmids):
    """Return {pmid: {'title':..., 'year':...}} — order-independent."""
    if not pmids:
        return {}

    handle = Entrez.efetch(db="pubmed", id=",".join(pmids),
                           rettype="medline", retmode="text")
    out = {}
    for art in Medline.parse(handle):
        pmid = (art.get("PMID") or "").strip()
        if not pmid:
            continue
        out[pmid] = {"title": art.get("TI", "").strip(),
                     "year":  extract_year(art.get("DP", ""))}
    handle.close()
    return out


# ## 6 · Main screen

# In[ ]:


def run_screen(genes, alias_dict, retmax=RETMAX, sort=SORT, verbose=True):
    disease_block = build_disease_block()
    cache, rows = {}, []     

    for gene in genes:
        gene = str(gene).strip()

        if gene in cache:                                
            rows.append({**cache[gene], "gene_name": gene})
            continue

        aliases = alias_dict.get(gene, [gene])
        row = {
            "gene_name": gene,
            "n_aliases": len(aliases),
            "aliases_used": "; ".join(aliases),
            "paper_count": None,
            "top_pmids": "", "top_titles": "", "publication_years": "",
            "query": "", "query_translation": "", "api_flags": "",
            "status": "not_run",
        }

        try:
            query = build_query(aliases, disease_block)
            row["query"] = query

            res = esearch_checked(query, retmax=retmax, sort=sort)
            row["paper_count"]       = res["count"]
            row["top_pmids"]         = "; ".join(res["pmids"])
            row["query_translation"] = res["query_translation"]
            row["api_flags"]         = res["flags"]

            details = efetch_by_pmid(res["pmids"])          
            row["top_titles"] = " || ".join(
                details.get(p, {}).get("title", "<not returned>") for p in res["pmids"]
            )
            row["publication_years"] = "; ".join(
                details.get(p, {}).get("year", "") for p in res["pmids"]
            )

            n_missing = len([p for p in res["pmids"] if p not in details])
            row["status"] = "ok" if res["clean"] else "SUSPECT_QUERY"
            if n_missing:
                row["status"] += f" (+{n_missing} PMIDs not returned by efetch)"

        except Exception as e:
            row["status"] = f"FAILED: {type(e).__name__}: {e}"
            if verbose:
                print(f"  !! {gene}: {row['status']}", file=sys.stderr)

        cache[gene] = row
        rows.append(row)



    df = pd.DataFrame(rows)

    df = df.sort_values(
        by=["paper_count"], ascending=True, na_position="first"
    ).reset_index(drop=True)

    return df


# ## 7 · Verification

# In[ ]:


def _check(name, condition, detail=""):
    mark = "PASS" if condition else "FAIL"
    print(f"  [{mark}] {name}" + (f"  -- {detail}" if detail and not condition else ""))
    return bool(condition)


def selftest_offline():
    print("\n=== GROUP A: OFFLINE SELF-TESTS ===")
    ok = []
    approved = {"ERBB2", "NEU4", "TP53", "LOX", "MLH1"}

    keep, reason = sanitize_alias("p185(erbB2)", "ERBB2", approved)
    ok.append(_check("A1  'p185(erbB2)' is rejected (unsafe parens)", keep is False, reason))

    keep, _ = sanitize_alias("LOX", "LOX", approved)
    ok.append(_check("A2  approved symbol 'LOX' survives despite blocklist", keep is True))


    keep, _ = sanitize_alias("LOX", "SOMEGENE", approved)
    ok.append(_check("A3  'LOX' rejected as a foreign alias", keep is False))

    keep, _ = sanitize_alias("NEU", "ERBB2", approved)
    ok.append(_check("A4  3-letter alias 'NEU' rejected", keep is False))


    keep, _ = sanitize_alias("CD340", "ERBB2", approved)
    ok.append(_check("A5  'CD340' kept", keep is True))


    erbb2_raw = ['CD340','ERBB2','HER-2','HER2','MLN-19','NEU',
                 'c-ERB-2','c-ERB2','p185(erbB2)']
    old_q = build_query_ORIGINAL(erbb2_raw)
    ok.append(_check(
        "A6  ORIGINAL query is paren-BALANCED (so balance-checking is useless)",
        old_q.count("(") == old_q.count(")"),
        f"{old_q.count('(')} open / {old_q.count(')')} close",
    ))
    print(f"       offending term: ...{old_q[old_q.find('p185')-10:old_q.find('p185')+28]}...")


    try:
        assert_query_wellformed(old_q)
        ok.append(_check("A7  ORIGINAL query rejected by validator", False, "it passed!"))
    except ValueError as e:
        ok.append(_check("A7  ORIGINAL query rejected by validator", True))
        for line in str(e).splitlines()[0].split("; ")[:3]:
            print(f"       reason: {line[:110]}")


    safe = [a for a in erbb2_raw if sanitize_alias(a, "ERBB2", approved)[0]]
    new_q = build_query(safe)
    try:
        assert_query_wellformed(new_q)
        ok.append(_check("A8  CORRECTED query passes validator", True))
    except ValueError as e:
        ok.append(_check("A8  CORRECTED query passes validator", False, str(e)))
    print(f"       kept aliases: {safe}")


    faults = {
        '(A[tiab] AND (B[tiab])':        "unbalanced parens",
        'A[tiab] AND NOT review[pt]':    "adjacent operators",
        '"A[tiab] AND B[tiab]':          "odd quotes",
        '() AND B[tiab]':                "empty group",
        '"p185(erbB2)"[tiab] OR "X"[tiab]': "balanced-but-regrouping quoted phrase",
        ')A[tiab] AND B[tiab](':         "close before open",
    }
    caught = 0
    for bad, label in faults.items():
        try:
            assert_query_wellformed(bad)
        except ValueError:
            caught += 1
    ok.append(_check(f"A9  validator caught {caught}/{len(faults)} planted faults",
                     caught == len(faults)))

    cases = [("2016 Jan 15", "2016"), ("Winter 2015", "2015"),
             ("2016 Nov-Dec", "2016"), ("", ""), ("2020", "2020")]
    good = all(extract_year(a) == b for a, b in cases)
    ok.append(_check("A10 extract_year handles awkward DP strings", good,
                     str([(a, extract_year(a)) for a, b in cases if extract_year(a) != b])))

    fixture = textwrap.dedent("""\
        PMID- 33333333
        TI  - Third paper returned first, out of order
        DP  - 2019 Nov-Dec

        PMID- 11111111
        TI  - A very long title that wraps across
              two physical lines in MEDLINE format
        DP  - 2016 Jan 15

        PMID- 22222222
        DP  - 2020 Mar
        """)
    recs = {}
    for art in Medline.parse(io.StringIO(fixture)):
        pmid = (art.get("PMID") or "").strip()
        if pmid:
            recs[pmid] = {"title": art.get("TI", "").strip(),
                          "year": extract_year(art.get("DP", ""))}
    requested = ["11111111", "22222222", "33333333", "99999999"]  # last one absent
    titles = [recs.get(p, {}).get("title", "<not returned>") for p in requested]
    ok.append(_check("A11 titles map to the RIGHT pmid despite reordering",
                     titles[0].startswith("A very long title")))
    ok.append(_check("A11b record with no TI does not shift the others",
                     titles[1] == "" and titles[2].startswith("Third paper")))
    ok.append(_check("A11c missing PMID is flagged, not silently dropped",
                     titles[3] == "<not returned>"))

    ok.append(_check("A12 pipe delimiter detected",
                     detect_multivalue_delimiter(pd.Series(["A|B|C", "D|E"])) == "|"))
    ok.append(_check("A12b comma delimiter detected",
                     detect_multivalue_delimiter(pd.Series(["A, B", "C, D"])) == ","))

    print(f"\n  OFFLINE RESULT: {sum(ok)}/{len(ok)} checks passed")
    return all(ok)


def selftest_online_controls():
    """Group B1 — positive / negative / structural controls against live PubMed.
    Requires network access to eutils.ncbi.nlm.nih.gov."""
    print("\n=== GROUP B1: LIVE PUBMED CONTROLS ===")
    disease = build_disease_block()
    ok = []

    r = esearch_checked(build_query(["TP53"], disease))
    ok.append(_check(f"B1  positive control TP53 count={r['count']} (expect >500)",
                     r["count"] > 500))
    print(f"       QueryTranslation: {r['query_translation'][:160]}")

    r = esearch_checked(build_query(["ZZQQXNOTAGENE"], disease))
    ok.append(_check(f"B2  negative control count={r['count']} (expect 0)", r["count"] == 0))
    ok.append(_check("B3  negative control raises a PubMed warning", bool(r["flags"]),
                     "no flags returned"))
    print(f"       flags: {r['flags'][:200]}")

    a = esearch_checked(build_query(["ERBB2"], disease))["count"]
    b = esearch_checked(build_query(["ERBB2", "CD340"], disease))["count"]
    ok.append(_check(f"B4  OR-ing an alias never lowers count ({a} -> {b})", b >= a))

    g = esearch_checked('"ERBB2"[Title/Abstract]')["count"]
    ok.append(_check(f"B5  gene-alone ({g}) > gene AND disease ({a})", g > a))

  
    c1 = esearch_checked(build_query(["ERBB2"], disease), retmax=1)["count"]
    ok.append(_check(f"B6  Count independent of retmax ({a} vs {c1})", a == c1))

    print("\n  --- OLD vs NEW on ERBB2 (the diagnostic comparison) ---")
    erbb2_raw = ['CD340','ERBB2','HER-2','HER2','MLN-19','NEU',
                 'c-ERB-2','c-ERB2','p185(erbB2)']
    try:
        h = Entrez.esearch(db="pubmed", term=build_query_ORIGINAL(erbb2_raw),
                           retmax=1, sort=SORT)
        old = Entrez.read(h); h.close()
        print(f"       ORIGINAL query Count = {old['Count']}")
        print(f"       ORIGINAL translation = {old.get('QueryTranslation','')[:200]}")
        print(f"       ORIGINAL warnings    = {dict(old.get('WarningList', {}))}")
    except Exception as e:
        print(f"       ORIGINAL query raised: {type(e).__name__}: {e}")
    print(f"       CORRECTED Count      = {a}")

    print(f"\n  LIVE CONTROL RESULT: {sum(ok)}/{len(ok)} checks passed")
    return all(ok)


def crosscheck_europepmc(gene, aliases=None, timeout=30):
    import requests
    aliases = aliases or [gene]
    gene_block = " OR ".join(f'TITLE_ABS:"{a}"' for a in aliases)
    disease_block = " OR ".join(
        f'TITLE_ABS:"{t}"' for t in
        ["esophageal cancer", "oesophageal cancer", "esophageal carcinoma",
         "esophageal squamous cell carcinoma", "esophageal adenocarcinoma"]
    )
    q = f"({gene_block}) AND ({disease_block})"
    r = requests.get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        params={"query": q, "format": "json", "pageSize": 1, "synonym": "false"},
        timeout=timeout,
    )
    r.raise_for_status()
    return {"gene": gene, "europepmc_hitCount": r.json().get("hitCount"), "query": q}


def crosscheck_pubtator3(gene, timeout=30):

    import requests
    r = requests.get(
        "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/",
        params={"text": f"@GENE_{gene} AND @DISEASE_Esophageal_Neoplasms"},
        timeout=timeout,
    )
    r.raise_for_status()
    j = r.json()
    return {"gene": gene, "pubtator3_count": j.get("count"), "raw_keys": list(j)[:8]}


def audit_report(df, alias_audit, out_prefix="litscreen"):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    print("\n=== AUDIT SUMMARY ===")
    print(f"  run timestamp (UTC) : {stamp}")
    print(f"  genes screened      : {len(df)}")
    print(f"  status == ok        : {(df['status'] == 'ok').sum()}")
    print(f"  SUSPECT_QUERY       : {df['status'].str.startswith('SUSPECT').sum()}")
    print(f"  FAILED              : {df['status'].str.startswith('FAILED').sum()}")
    print(f"  zero-hit genes      : {(df['paper_count'] == 0).sum()}")

    bad = df[df["status"] != "ok"]
    if len(bad):
        print("\n  !! Genes needing manual review:")
        for _, r in bad.iterrows():
            print(f"     {r['gene_name']:<10} {r['status']}")
            if r["api_flags"]:
                print(f"                {r['api_flags'][:150]}")

    df.to_csv(f"{out_prefix}_results.csv", index=False)
    alias_audit.to_csv(f"{out_prefix}_alias_audit.csv", index=False)
    with open(f"{out_prefix}_runmeta.json", "w") as fh:
        json.dump({"timestamp_utc": stamp, "retmax": RETMAX, "sort": SORT,
                   "exclude_reviews": EXCLUDE_REVIEWS,
                   "alias_min_len": ALIAS_MIN_LEN,
                   "disease_block": build_disease_block(),
                   "api_key_used": bool(Entrez.api_key),
                   "biopython": __import__("Bio").__version__}, fh, indent=2)
    print(f"\n  wrote {out_prefix}_results.csv, {out_prefix}_alias_audit.csv, "
          f"{out_prefix}_runmeta.json")


# ---
# # 8 · VERIFICATION

# In[ ]:


selftest_offline()


# In[ ]:


hgnc, delim = load_hgnc("custom.tsv")
alias_dict, alias_audit = build_alias_dict(hgnc, delim)

print(f"\nBuilt alias sets for {len(alias_dict):,} approved genes")
print(f"Aliases kept   : {alias_audit.kept.sum():,}")
print(f"Aliases dropped: {(~alias_audit.kept).sum():,}\n")
print(alias_audit[~alias_audit.kept].reason.value_counts().to_string())

print("\nERBB2 before/after:")
print("  BEFORE (yours):", ['CD340','ERBB2','HER-2','HER2','MLN-19','NEU',
                            'c-ERB-2','c-ERB2','p185(erbB2)'])
print("  AFTER  (fixed):", alias_dict.get("ERBB2"))


# In[ ]:


selftest_online_controls()


# In[ ]:


SUSPECT = ["TYMS", "RAD51", "PLK1", "MLH1", "PCNA", "HRAS"]

rows = []
for g in SUSPECT:
    aliases = alias_dict.get(g, [g])
    raw = [g] + split_multivalue(
        hgnc.loc[hgnc["Approved symbol"] == g, "Alias symbols"].squeeze()
        if (hgnc["Approved symbol"] == g).any() else "", delim)

    try:
        h = Entrez.esearch(db="pubmed", term=build_query_ORIGINAL(raw),
                           retmax=1, sort=SORT)
        old = Entrez.read(h); h.close()
        old_count, old_warn = int(old["Count"]), dict(old.get("WarningList", {}))
    except Exception as e:
        old_count, old_warn = None, f"{type(e).__name__}: {e}"

    new = esearch_checked(build_query(aliases))

    rows.append({"gene": g, "OLD_count": old_count, "NEW_count": new["count"],
                 "old_warnings": str(old_warn)[:90], "new_flags": new["flags"][:90]})
    print(f"{g:<7} OLD={old_count!s:<7} NEW={new['count']:<7}")
    print(f"        NEW translation: {new['query_translation'][:150]}")

import pandas as pd
pd.DataFrame(rows)


# In[ ]:


novelnew = pd.read_csv("try_file.csv")
novelnew = novelnew.dropna(subset=["gene_name"]).copy()

results = run_screen(novelnew["gene_name"].astype(str).str.strip().tolist(),
                     alias_dict)

merged = novelnew.merge(
    results.drop_duplicates("gene_name"),
    on="gene_name", how="left", suffixes=("", "_lit"))

merged = merged.sort_values(
    by=["paper_count", "class_1_mean"],
    ascending=[True, False], na_position="first").reset_index(drop=True)

audit_report(merged, alias_audit, out_prefix="ESCA_litscreen")
merged.head(25)


# In[ ]:


check = results.head(10)["gene_name"].tolist()

cross = []
for g in check:
    row = {"gene": g,
           "pubmed": int(results.loc[results.gene_name == g, "paper_count"].iloc[0])}
    try:
        row["europepmc"] = crosscheck_europepmc(g, alias_dict.get(g, [g]))["europepmc_hitCount"]
    except Exception as e:
        row["europepmc"] = f"ERR {type(e).__name__}"
    try:
        row["pubtator3"] = crosscheck_pubtator3(g)["pubtator3_count"]
    except Exception as e:
        row["pubtator3"] = f"ERR {type(e).__name__}"
    cross.append(row); print(row)

cross_df = pd.DataFrame(cross)
cross_df["ratio_epmc"] = pd.to_numeric(cross_df.europepmc, errors="coerce") / \
                         cross_df.pubmed.replace(0, pd.NA)
print("\nFlag any gene where ratio_epmc > 10 or < 0.1 — its query is suspect.")
cross_df

