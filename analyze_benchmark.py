#!/usr/bin/env python3
"""Tabela comparativa de benchmarks SpMV entre backends (matrizes em comum)."""
import argparse
import io
import sys
import pandas as pd


def load(path):
    with open(path, encoding="utf-8-sig") as f:
        lines = [ln.rstrip().rstrip(",").rstrip() for ln in f if ln.strip()]
    raw = pd.read_csv(io.StringIO("\n".join(lines)),
                      header=None, dtype=str, skipinitialspace=True)

    hdr = [str(c).strip() for c in raw.iloc[0]]
    data = raw.iloc[1:].reset_index(drop=True)
    while len(hdr) < data.shape[1]:
        hdr.append("Max_err_rel" if len(hdr) == data.shape[1] - 1
                   else f"Unnamed_{len(hdr)}")
    data.columns = hdr[:data.shape[1]]

    for c in data.columns:
        if str(c).lower() in ("problema", "problem", "name"):
            data = data.rename(columns={c: "Problema"})
            break

    data["Problema"] = data["Problema"].astype(str).str.strip()
    for c in data.columns:
        if c != "Problema":
            data[c] = pd.to_numeric(data[c], errors="coerce")
    return data


# (nome_coluna, header, formatador) — campos exibidos por backend
COLS = [
    ("nnz",                "nnz",      lambda v: f"{int(v)}"),
    ("Dim",                "Dim",      lambda v: f"{int(v)}"),
    ("niter",              "niter",    lambda v: f"{int(v)}"),
    ("t_GPU_copy_ms",      "copy",     lambda v: f"{v:.3e}"),
    ("t_GPU_prep_ms",      "prep",     lambda v: f"{v:.3e}"),
    ("t_GPU_SpMV_avg_ms",  "gpu_avg",  lambda v: f"{v:.3e}"),
    ("t_CPU_SpMV_avg_ms",  "cpu_avg",  lambda v: f"{v:.3e}"),
    ("Speedup",            "sp",       lambda v: f"{v:.3f}"),
    ("Speedup_with_copy",  "sp+copy",  lambda v: f"{v:.3f}"),
    ("Max_err_rel",        "err",      lambda v: f"{v:.2e}"),
]

W_FIELD = 11
W_NAME  = 22
W_NNZ   = 10
W_DIM   = 8
W_BACK  = len(COLS) * W_FIELD


def auto_label(path):
    """Deriva rotulo curto do nome do arquivo."""
    name = path.rsplit("/", 1)[-1].replace(".csv", "")
    low = name.lower()
    arch = "NVIDIA" if "nvidia" in low else "AMD" if "amd" in low else "?"
    kind = "legacy" if "legacy" in low else "generic"
    return f"{arch} {kind}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csvs", nargs="+", help="arquivos CSV")
    ap.add_argument("--labels", nargs="+", default=None,
                    help="rotulos manuais (um por CSV)")
    args = ap.parse_args()

    labels = args.labels or [auto_label(p) for p in args.csvs]
    if len(labels) != len(args.csvs):
        print("--labels deve ter o mesmo numero de CSVs", file=sys.stderr)
        sys.exit(1)

    dfs = [load(p) for p in args.csvs]
    common = set.intersection(*[set(d["Problema"]) for d in dfs])
    if not common:
        print("Sem matrizes em comum.")
        return

    idx = {l: d.set_index("Problema") for l, d in zip(labels, dfs)}
    names = sorted(common, key=lambda n: idx[labels[0]].loc[n, "nnz"])

    # --- cabecalho linha 1: rotulos centralizados sobre o bloco ---
    h1 = f"{'Matriz':<{W_NAME}} {'nnz':>{W_NNZ}} {'Dim':>{W_DIM}}"
    for l in labels:
        h1 += "|" + l.center(W_BACK)
    h1 += "|"

    # --- cabecalho linha 2: nomes das metricas ---
    h2 = " " * (W_NAME + 1 + W_NNZ + 1 + W_DIM)
    for _ in labels:
        h2 += "|"
        for _, short, _ in COLS:
            h2 += f"{short:>{W_FIELD}}"
    h2 += "|"

    sep = "-" * max(len(h1), len(h2))
    print(h1)
    print(sep)
    print(h2)
    print(sep)

    # --- dados ---
    for name in names:
        r0 = idx[labels[0]].loc[name]
        row = (f"{name:<{W_NAME}} {int(r0['nnz']):>{W_NNZ}} "
               f"{int(r0['Dim']):>{W_DIM}}")
        for l in labels:
            row += "|"
            r = idx[l].loc[name]
            for col, _, fmt in COLS:
                if col in r.index and pd.notna(r[col]):
                    s = fmt(r[col])
                else:
                    s = "—"
                row += f"{s:>{W_FIELD}}"
        row += "|"
        print(row)

    print()
    print(f"Matrizes em comum: {len(common)}")


if __name__ == "__main__":
    main()