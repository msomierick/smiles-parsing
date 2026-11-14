from smiles_parsing.__main__ import main


def test_prints_one_token_per_line(capsys):
    assert main("CO") == 0
    assert capsys.readouterr().out.splitlines() == [
        "ELEMENT        'C'      (0, 1)",
        "ELEMENT        'O'      (1, 2)",
    ]


def test_prints_errors_to_stderr(capsys):
    assert main("CXC") == 1
    assert capsys.readouterr().err == "unknown atom symbol 'X'\nCXC\n ^\n"
