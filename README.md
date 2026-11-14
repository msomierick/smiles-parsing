# smiles-parsing

Companion code for the [Parsing SMILES](https://labofcoding.com/posts/parsing-smiles-from-smiles-strings-to-molecular-graphs) blog series, which builds a SMILES parser in Python by hand, one stage per post. Each post ends at a git tag, so you can check out the code as it stands at the end of that post:

- [Lexical analysis](https://labofcoding.com/posts/parsing-smiles-lexical-analysis): tag `part-2-lexer`

## Setup

You need [git](https://git-scm.com/) and [uv](https://docs.astral.sh/uv/). Clone the repository and check out the tag for the post you are reading:

```bash
git clone https://github.com/msomierick/smiles-parsing.git
cd smiles-parsing
git checkout part-2-lexer
```

Then, from the repository root:

```bash
uv run python -m smiles_parsing "[13CH4]"   # print the tokens of a SMILES string
uv run pytest                              # run the tests
uv run ruff check                          # lint
```

uv creates the virtual environment and installs the dependencies the first time you run it.
