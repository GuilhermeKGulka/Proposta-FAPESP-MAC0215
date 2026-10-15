#!/usr/bin/env python3
"""Gera a lista de matrizes (matrices.txt) a partir de filtros, como o ssgui.

Exemplo (o que aparece na janela do ssgui: quadrada, SPD, real, simetria numerica 1.0):
    python select_matrices.py --square --spd yes --dtype real --nsym-min 1 --nsym-max 1 -o matrices.txt

Depois: python download_matrices.py matrices.txt   (opcionalmente com --max-nnz)

Os filtros sao aplicados aqui, sobre o indice completo da colecao (o mesmo que o
add_matrix_kind.py usa). 'python select_matrices.py --show-fields' mostra os
atributos disponiveis em cada matriz.
"""
import argparse
import sys

import ssgetpy

# Em algumas versoes do ssgetpy, limit=0 devolve ZERO resultados; usa-se um limite alto.
LIMIT_ALL = 100000


def yesno(v):
    return {"yes": True, "no": False}[v]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--out", default="matrices.txt", help="arquivo de saida (padrao: matrices.txt)")
    ap.add_argument("--rows-min", type=int); ap.add_argument("--rows-max", type=int)
    ap.add_argument("--cols-min", type=int); ap.add_argument("--cols-max", type=int)
    ap.add_argument("--nnz-min", type=int);  ap.add_argument("--nnz-max", type=int)
    ap.add_argument("--square", action="store_true", help="so matrizes quadradas")
    ap.add_argument("--spd", choices=["yes", "no"], help="positiva definida?")
    ap.add_argument("--discretization", choices=["yes", "no"], help="discretizacao 2D/3D?")
    ap.add_argument("--dtype", choices=["real", "complex", "binary"])
    ap.add_argument("--psym-min", type=float); ap.add_argument("--psym-max", type=float)
    ap.add_argument("--nsym-min", type=float); ap.add_argument("--nsym-max", type=float)
    ap.add_argument("--group", action="append", help="grupo (pode repetir)")
    ap.add_argument("--kind", action="append", help="trecho do 'kind' (pode repetir; vale qualquer um)")
    ap.add_argument("--show-fields", action="store_true", help="mostra os atributos de uma matriz e sai")
    args = ap.parse_args()

    print("Carregando indice da SuiteSparse (a 1a vez baixa o indice)...")
    allm = list(ssgetpy.search(limit=LIMIT_ALL))
    if not allm:
        sys.exit("O indice da SuiteSparse veio vazio. Verifique a rede e o cache do ssgetpy "
                 "(apague o cache para baixar o indice de novo). Teste:\n"
                 "  python -c \"import ssgetpy; print(len(ssgetpy.search(limit=5)), "
                 "len(ssgetpy.search(name='bcsstk01')))\"")
    if args.show_fields:
        print(vars(allm[0]) if hasattr(allm[0], "__dict__") else dir(allm[0]))
        return 0

    def get(m, attr):
        try:
            return getattr(m, attr)
        except AttributeError:
            sys.exit(f"Atributo '{attr}' nao existe nas matrizes desta versao do ssgetpy; "
                     f"veja --show-fields.")

    # lista de (descricao, predicado), aplicada em sequencia para mostrar o efeito de cada filtro
    checks = []
    def add(desc, pred):
        checks.append((desc, pred))

    if args.rows_min is not None: add(f"rows >= {args.rows_min}", lambda m: get(m, "rows") >= args.rows_min)
    if args.rows_max is not None: add(f"rows <= {args.rows_max}", lambda m: get(m, "rows") <= args.rows_max)
    if args.cols_min is not None: add(f"cols >= {args.cols_min}", lambda m: get(m, "cols") >= args.cols_min)
    if args.cols_max is not None: add(f"cols <= {args.cols_max}", lambda m: get(m, "cols") <= args.cols_max)
    if args.nnz_min is not None:  add(f"nnz >= {args.nnz_min}", lambda m: get(m, "nnz") >= args.nnz_min)
    if args.nnz_max is not None:  add(f"nnz <= {args.nnz_max}", lambda m: get(m, "nnz") <= args.nnz_max)
    if args.square:               add("quadrada", lambda m: get(m, "rows") == get(m, "cols"))
    if args.spd is not None:      add(f"spd = {args.spd}", lambda m: bool(get(m, "isspd")) == yesno(args.spd))
    if args.discretization is not None:
        add(f"2D/3D = {args.discretization}", lambda m: bool(get(m, "is2d3d")) == yesno(args.discretization))
    if args.dtype is not None:    add(f"dtype = {args.dtype}", lambda m: str(get(m, "dtype")).lower() == args.dtype)
    if args.psym_min is not None: add(f"psym >= {args.psym_min}", lambda m: get(m, "psym") >= args.psym_min)
    if args.psym_max is not None: add(f"psym <= {args.psym_max}", lambda m: get(m, "psym") <= args.psym_max)
    if args.nsym_min is not None: add(f"nsym >= {args.nsym_min}", lambda m: get(m, "nsym") >= args.nsym_min)
    if args.nsym_max is not None: add(f"nsym <= {args.nsym_max}", lambda m: get(m, "nsym") <= args.nsym_max)
    if args.group:                add(f"grupo em {args.group}", lambda m: get(m, "group") in args.group)
    if args.kind:                 add(f"kind contem {args.kind}",
                                      lambda m: any(k.lower() in str(get(m, "kind")).lower() for k in args.kind))

    current = list(allm)
    print(f"{len(current):6d}  (indice completo)")
    for desc, pred in checks:
        current = [m for m in current if pred(m)]
        print(f"{len(current):6d}  apos filtro: {desc}")

    if not current and checks:
        print("\nNenhuma matriz passou. Valores encontrados no indice (confira a representacao):")
        for attr in ("dtype", "isspd", "is2d3d"):
            if hasattr(allm[0], attr):
                print(f"  {attr}: {sorted({repr(getattr(m, attr)) for m in allm})[:10]}")
        for attr in ("psym", "nsym"):
            if hasattr(allm[0], attr):
                vals = [getattr(m, attr) for m in allm]
                print(f"  {attr}: min={min(vals)!r} max={max(vals)!r}")

    chosen = sorted(current, key=lambda m: (m.group.lower(), m.name.lower()))
    with open(args.out, "w", encoding="utf-8") as f:
        f.write("# Gerado por select_matrices.py: " + " ".join(sys.argv[1:]) + "\n")
        f.write(f"# {len(chosen)} matrizes\n")
        for m in chosen:
            f.write(f"{m.group}/{m.name}  # {m.rows}x{m.cols} nnz={m.nnz} {m.kind}\n")
    print(f"{len(chosen)} matrizes selecionadas -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())