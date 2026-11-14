import sys

from smiles_parsing.lexer import LexicalError, lex


def main(smiles: str) -> int:
    try:
        tokens = lex(smiles)
    except LexicalError as error:
        print(error, file=sys.stderr)
        return 1
    for token in tokens:
        print(f"{token.type.name:<14} {token.lexeme!r:<8} {token.span}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python -m smiles_parsing SMILES")
    sys.exit(main(sys.argv[1]))
