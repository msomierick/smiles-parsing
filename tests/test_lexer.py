import pytest

from smiles_parsing.lexer import ELEMENTS, LexicalError, lex
from smiles_parsing.tokens import TokenType as T


def lexemes(source: str) -> list[str]:
    return [token.lexeme for token in lex(source)]


def types(source: str) -> list[T]:
    return [token.type for token in lex(source)]


# Inputs the lexer can tokenize. Some, like C(C, are not valid SMILES: later
# stages reject them.
TOKENIZABLE = [
    ("CCO", ["C", "C", "O"]),
    ("c1ccccc1", ["c", "1", "c", "c", "c", "c", "c", "1"]),
    ("C12CC12", ["C", "1", "2", "C", "C", "1", "2"]),
    ("CC(=O)O", ["C", "C", "(", "=", "O", ")", "O"]),
    ("F/C=C/F", ["F", "/", "C", "=", "C", "/", "F"]),
    ("C#N", ["C", "#", "N"]),
    ("C.C", ["C", ".", "C"]),
    ("*C", ["*", "C"]),
    ("", []),
    ("C(C", ["C", "(", "C"]),
    ("[C@TB99]", ["[", "C", "@TB99", "]"]),
    ("[C+H]", ["[", "C", "+", "H", "]"]),
    # Maximum munch
    ("ClCCl", ["Cl", "C", "Cl"]),
    ("BrCB", ["Br", "C", "B"]),
    ("C%12CC%12", ["C", "%12", "C", "C", "%12"]),
    ("C%123", ["C", "%12", "3"]),
    # Only organic subset symbols outside brackets: Cs is C then aromatic s
    ("Cs", ["C", "s"]),
    ("Sc", ["S", "c"]),
    ("Co", ["C", "o"]),
    # Bracket atoms
    ("[13CH4]", ["[", "13", "C", "H4", "]"]),
    ("[2H]", ["[", "2", "H", "]"]),
    ("[O-]C(=O)C", ["[", "O", "-", "]", "C", "(", "=", "O", ")", "C"]),
    ("[Fe+2]", ["[", "Fe", "+2", "]"]),
    ("[Ti-12]", ["[", "Ti", "-12", "]"]),
    ("[Cu++]", ["[", "Cu", "++", "]"]),
    ("[CH3:1]", ["[", "C", "H3", ":1", "]"]),
    ("[NH4+:12]", ["[", "N", "H4", "+", ":12", "]"]),
    ("[HH]", ["[", "H", "H", "]"]),
    ("[Na+].[Cl-]", ["[", "Na", "+", "]", ".", "[", "Cl", "-", "]"]),
    ("[Sc]", ["[", "Sc", "]"]),
    ("[Hg]", ["[", "Hg", "]"]),
    ("[se]", ["[", "se", "]"]),
    ("c1cc[nH]c1", ["c", "1", "c", "c", "[", "n", "H", "]", "c", "1"]),
    ("[*]", ["[", "*", "]"]),
    (
        "N[C@@H](C)C(=O)O",
        ["N", "[", "C", "@@", "H", "]", "(", "C", ")", "C", "(", "=", "O", ")", "O"],
    ),
    ("[C@TH1]", ["[", "C", "@TH1", "]"]),
    ("[C@AL2]", ["[", "C", "@AL2", "]"]),
    ("[Pt@SP3]", ["[", "Pt", "@SP3", "]"]),
    ("[As@TB7]", ["[", "As", "@TB7", "]"]),
    ("[Fe@OH12]", ["[", "Fe", "@OH12", "]"]),
]


@pytest.mark.parametrize(("source", "expected"), TOKENIZABLE)
def test_lexemes(source, expected):
    assert lexemes(source) == expected


@pytest.mark.parametrize("source", [source for source, _ in TOKENIZABLE])
def test_spans_match_lexemes(source):
    for token in lex(source):
        start, end = token.span
        assert source[start:end] == token.lexeme


@pytest.mark.parametrize("source", [source for source, _ in TOKENIZABLE])
def test_lexemes_rebuild_the_source(source):
    assert "".join(token.lexeme for token in lex(source)) == source


def test_token_types():
    assert types("c1C(=O)*.") == [
        T.AROMATIC,
        T.RING_NUMBER,
        T.ELEMENT,
        T.BRANCH_START,
        T.BOND,
        T.ELEMENT,
        T.BRANCH_END,
        T.WILDCARD,
        T.DOT,
    ]


def test_bracket_token_types():
    assert types("[13CH4+:2]") == [
        T.BRACKET_START,
        T.ISOTOPE,
        T.ELEMENT,
        T.HYDROGEN_COUNT,
        T.CHARGE,
        T.ATOM_CLASS,
        T.BRACKET_END,
    ]
    assert types("[C@@H-]") == [
        T.BRACKET_START,
        T.ELEMENT,
        T.CHIRALITY,
        T.HYDROGEN_COUNT,
        T.CHARGE,
        T.BRACKET_END,
    ]


def test_h_is_an_element_before_the_symbol_and_a_count_after_it():
    assert types("[2H]") == [T.BRACKET_START, T.ISOTOPE, T.ELEMENT, T.BRACKET_END]
    assert types("[CH]") == [
        T.BRACKET_START,
        T.ELEMENT,
        T.HYDROGEN_COUNT,
        T.BRACKET_END,
    ]


def test_every_element_is_known():
    assert len(ELEMENTS) == 118


@pytest.mark.parametrize(
    ("source", "span"),
    [
        ("CXC", (1, 2)),
        ("C?", (1, 2)),
        ("Cé", (1, 2)),
        ("[é]", (1, 2)),
        ("H", (0, 1)),
        ("Na", (1, 2)),
        ("c1ccsec1", (4, 6)),
        ("c1ccasc1", (4, 6)),
        ("C%", (1, 2)),
        ("C%1", (1, 3)),
        ("C%1C", (1, 3)),
        ("C]", (1, 2)),
        ("[C", (0, 2)),
        ("C[", (1, 2)),
        ("[Q]", (1, 2)),
        ("[C[C]]", (2, 3)),
        ("[C@TH]", (2, 5)),
        ("[C@TH3]", (2, 5)),
        ("[C@SP4]", (2, 5)),
        ("[C@TB]", (2, 5)),
        ("[C H4]", (2, 3)),
        ("[C1]", (2, 3)),
        ("[CC]", (2, 3)),
        ("[C*]", (2, 3)),
        ("[CHe]", (3, 4)),
        ("[CH12]", (2, 5)),
        ("[C:]", (2, 3)),
        ("[Fe@OH123]", (8, 9)),
        ("[Fe+123]", (6, 7)),
        ("[Fe+++]", (3, 6)),
        ("[Fe+-]", (3, 5)),
        ("[Fe+2-]", (3, 6)),
    ],
)
def test_errors(source, span):
    with pytest.raises(LexicalError) as error:
        lex(source)
    assert error.value.span == span


def test_error_points_at_the_problem():
    with pytest.raises(LexicalError) as error:
        lex("CXC")
    assert str(error.value) == "unknown atom symbol 'X'\nCXC\n ^"


def test_error_underlines_the_whole_span():
    with pytest.raises(LexicalError) as error:
        lex("C[CH4")
    assert str(error.value) == "bracket atom is not closed\nC[CH4\n ^^^^"


@pytest.mark.parametrize("terminator", [" ", "\t", "\n", "\r\n", "  "])
def test_stops_at_the_terminator(terminator):
    tokens = lex(f"CCO{terminator}ethanol")
    assert [token.lexeme for token in tokens] == ["C", "C", "O"]


def test_error_shows_only_the_first_line():
    with pytest.raises(LexicalError) as error:
        lex("[C\nX]")
    assert str(error.value) == "unexpected character '\\n' in bracket atom\n[C\n  ^"


def test_non_ascii_letters_are_unexpected_characters():
    with pytest.raises(LexicalError) as error:
        lex("Cé")
    assert error.value.message == "unexpected character 'é'"


def test_leading_whitespace_ends_an_empty_smiles():
    assert lex(" CCO") == []
