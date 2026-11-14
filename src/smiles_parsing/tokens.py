from dataclasses import dataclass
from enum import Enum, auto


class TokenType(Enum):
    # Atoms
    ELEMENT = auto()  # C, Cl, and any element symbol inside brackets
    AROMATIC = auto()  # c, n, and se or as inside brackets
    WILDCARD = auto()  # *

    # Bonds and structure
    BOND = auto()  # - = # $ : / \
    DOT = auto()  # .
    BRANCH_START = auto()  # (
    BRANCH_END = auto()  # )
    RING_NUMBER = auto()  # 1 or %12

    # Bracket atoms, such as [13CH4] or [NH4+:1]
    BRACKET_START = auto()  # [
    BRACKET_END = auto()  # ]
    ISOTOPE = auto()  # 13
    CHIRALITY = auto()  # @, @@, @TH1, @OH12
    HYDROGEN_COUNT = auto()  # H or H4
    CHARGE = auto()  # +, -, +2, -12, ++ or --
    ATOM_CLASS = auto()  # :1


@dataclass(frozen=True, slots=True)
class Token:
    type: TokenType
    lexeme: str
    span: tuple[int, int]  # half-open, so source[start:end] == lexeme
