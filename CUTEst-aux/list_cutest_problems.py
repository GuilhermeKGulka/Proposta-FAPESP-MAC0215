#!/usr/bin/env python3
"""
Gera uma lista de problemas CUTEst filtrados pela classificação oficial.

Os filtros aceitam as strings descritivas usadas pelo pycutest:

  objetivo:
    none | constant | linear | quadratic | sum of squares | other

  restrições:
    unconstrained | fixed | bound | adjacency | linear | quadratic | other

Uso:
    python list_cutest_problems.py
    python list_cutest_problems.py --constraints unconstrained --max-n 1000
    python list_cutest_problems.py --objective "sum of squares" --constraints unconstrained
    python list_cutest_problems.py --help
"""
import argparse
import sys

try:
    import pycutest
except ImportError:
    print("pycutest não encontrado. Instale com: pip install pycutest", file=sys.stderr)
    sys.exit(1)


OBJECTIVES = ["none", "constant", "linear", "quadratic", "sum of squares", "other"]
CONSTRAINTS = ["unconstrained", "fixed", "bound", "adjacency",
               "linear", "quadratic", "other"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Lista problemas CUTEst por classificação (padrão: unconstrained).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--objective", default=None, choices=OBJECTIVES,
                   help="Tipo de objetivo (padrão: sem filtro)")
    p.add_argument("--constraints", default="unconstrained", choices=CONSTRAINTS,
                   help="Tipo de restrição (padrão: unconstrained)")
    p.add_argument("--max-n", type=int, default=None,
                   help="Número máximo de variáveis (n)")
    p.add_argument("--max-m", type=int, default=None,
                   help="Número máximo de restrições gerais (m)")
    p.add_argument("-o", "--output", default="problems.txt",
                   help="Arquivo de saída (padrão: problems.txt)")
    return p.parse_args()


def main():
    args = parse_args()

    kwargs = {"constraints": args.constraints}
    if args.objective is not None:
        kwargs["objective"] = args.objective
    if args.max_n is not None:
        kwargs["n"] = (0, args.max_n)
    if args.max_m is not None:
        kwargs["m"] = (0, args.max_m)

    print(f"Filtrando CUTEst: {kwargs}...")
    names = sorted(pycutest.find_problems(**kwargs))

    with open(args.output, "w") as f:
        f.write("\n".join(names) + "\n")

    print(f"{len(names)} problemas encontrados; lista em {args.output}")


if __name__ == "__main__":
    main()