from smiles_parsing.tokens import Token, TokenType

DIGITS = "0123456789"
BONDS = "-=#$:/\\"
TERMINATORS = " \t\n\r"
ORGANIC_SUBSET = {
    "B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I",
    "b", "c", "n", "o", "p", "s",
}  # fmt: skip
ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
}  # fmt: skip
BRACKET_SYMBOLS = ELEMENTS | {"b", "c", "n", "o", "p", "s", "se", "as"}
# The grammar lists the digits allowed after @TH, @AL and @SP. @TB and @OH
# take one or two digits, and their ranges are checked in a later stage.
CHIRALITY_DIGITS = {"TH": "12", "AL": "12", "SP": "123"}


class LexicalError(Exception):
    def __init__(self, message: str, source: str, span: tuple[int, int]):
        super().__init__(message)
        self.message = message
        self.source = source
        self.span = span

    def __str__(self) -> str:
        # Errors happen at or before the first line break, so show only the first line
        line = self.source.splitlines()[0]
        start, end = self.span
        carets = " " * start + "^" * max(end - start, 1)
        return f"{self.message}\n{line}\n{carets}"


class Lexer:
    def __init__(self, source: str):
        self.source = source
        self.pos = 0
        self.tokens: list[Token] = []

    def lex(self) -> list[Token]:
        # The SMILES ends at the first terminator; what follows is not lexed
        while self.pos < len(self.source) and self._peek() not in TERMINATORS:
            self._lex_token()
        return self.tokens

    def _lex_token(self) -> None:
        char = self._peek()
        match char:
            case "(":
                self._add(TokenType.BRANCH_START, 1)
            case ")":
                self._add(TokenType.BRANCH_END, 1)
            case ".":
                self._add(TokenType.DOT, 1)
            case "*":
                self._add(TokenType.WILDCARD, 1)
            case "%":
                self._lex_percent_ring_number()
            case "[":
                self._lex_bracket_atom()
            case "]":
                # _lex_bracket_atom() reads every ']' that closes a bracket atom
                raise self._error("']' without a matching '['", 1)
            case _ if char in BONDS:
                self._add(TokenType.BOND, 1)
            case _ if char in DIGITS:
                self._add(TokenType.RING_NUMBER, 1)
            case _ if char.isascii() and char.isalpha():
                self._lex_atom(ORGANIC_SUBSET)
            case _:
                raise self._error(f"unexpected character {char!r}", 1)

    def _lex_bracket_atom(self) -> None:
        start = self.pos
        self._add(TokenType.BRACKET_START, 1)
        # Cases follow the grammar's order. Before the symbol, digits are an isotope
        # and letters start the symbol; after it, H is a hydrogen count
        seen_symbol = False
        while self.pos < len(self.source):
            char = self._peek()
            match char:
                case _ if char in DIGITS and not seen_symbol:
                    self._add(TokenType.ISOTOPE, self._count_digits(0))
                case "*" if not seen_symbol:
                    self._add(TokenType.WILDCARD, 1)
                    seen_symbol = True
                case _ if char.isascii() and char.isalpha() and not seen_symbol:
                    self._lex_atom(BRACKET_SYMBOLS)
                    seen_symbol = True
                case "@":
                    self._lex_chirality()
                case "H":
                    self._lex_hydrogen_count()
                case "+" | "-":
                    self._lex_charge()
                case ":":
                    self._lex_atom_class()
                case "]":
                    self._add(TokenType.BRACKET_END, 1)
                    return
                case _:
                    raise self._error(
                        f"unexpected character {char!r} in bracket atom", 1
                    )
        raise self._error("bracket atom is not closed", self.pos - start, start=start)

    def _lex_atom(self, symbols: set[str]) -> None:
        # Maximum munch: try the two-letter symbol first, so Cl is not C then l
        two_letters = self._peek() + self._peek(1)
        if two_letters in symbols:
            symbol = two_letters
        elif two_letters in ("se", "as"):
            raise self._error(f"{two_letters!r} is only allowed inside brackets", 2)
        elif self._peek() in symbols:
            symbol = self._peek()
        else:
            raise self._error(f"unknown atom symbol {self._peek()!r}", 1)
        token_type = TokenType.AROMATIC if symbol.islower() else TokenType.ELEMENT
        self._add(token_type, len(symbol))

    def _lex_hydrogen_count(self) -> None:
        digits = self._count_digits(1)
        if digits > 1:
            raise self._error("'H' can be followed by at most one digit", 1 + digits)
        self._add(TokenType.HYDROGEN_COUNT, 1 + digits)  # H or H4

    def _lex_charge(self) -> None:
        sign = self._peek()
        if self._peek(1) == sign:
            length = 2  # ++ or --
        else:
            length = 1 + min(self._count_digits(1), 2)  # +, -2
        if self._peek(length) in ("+", "-"):
            raise self._error("a charge cannot be followed by another sign", length + 1)
        self._add(TokenType.CHARGE, length)

    def _lex_atom_class(self) -> None:
        digits = self._count_digits(1)
        if digits == 0:
            raise self._error("':' must be followed by digits", 1)
        self._add(TokenType.ATOM_CLASS, 1 + digits)  # :1

    def _lex_percent_ring_number(self) -> None:
        digits = self._count_digits(1)
        if digits < 2:
            raise self._error("'%' must be followed by two digits", 1 + digits)
        self._add(TokenType.RING_NUMBER, 3)

    def _lex_chirality(self) -> None:
        if self._peek(1) == "@":
            self._add(TokenType.CHIRALITY, 2)  # @@
            return
        chirality_class = self._peek(1) + self._peek(2)
        if chirality_class in CHIRALITY_DIGITS:
            if self._peek(3) not in CHIRALITY_DIGITS[chirality_class]:
                allowed = " or ".join(CHIRALITY_DIGITS[chirality_class])
                raise self._error(
                    f"'@{chirality_class}' must be followed by {allowed}", 3
                )
            self._add(TokenType.CHIRALITY, 4)  # @TH1
        elif chirality_class in ("TB", "OH"):
            digits = min(self._count_digits(3), 2)
            if digits == 0:
                raise self._error(
                    f"'@{chirality_class}' must be followed by a digit", 3
                )
            self._add(TokenType.CHIRALITY, 3 + digits)  # @TB1 or @OH12
        else:
            self._add(TokenType.CHIRALITY, 1)  # @

    def _peek(self, offset: int = 0) -> str:
        """Return the character offset places ahead, or a NUL character past the end."""
        index = self.pos + offset
        return self.source[index] if index < len(self.source) else "\0"

    def _count_digits(self, offset: int) -> int:
        count = 0
        while self._peek(offset + count) in DIGITS:
            count += 1
        return count

    def _add(self, token_type: TokenType, length: int) -> None:
        end = self.pos + length
        self.tokens.append(
            Token(token_type, self.source[self.pos : end], (self.pos, end))
        )
        self.pos = end

    def _error(
        self, message: str, length: int, start: int | None = None
    ) -> LexicalError:
        if start is None:
            start = self.pos
        return LexicalError(message, self.source, (start, start + length))


def lex(source: str) -> list[Token]:
    return Lexer(source).lex()
