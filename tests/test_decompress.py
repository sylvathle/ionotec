import gzip
import pytest
from ionotec import decompress
from ionotec.download import extract


@pytest.mark.parametrize("module", [decompress, extract])
def test_extract_and_remove_gz(tmp_path, module):
    gz = tmp_path / "file.txt.gz"
    gz.write_bytes(gzip.compress(b"hello"))

    out = module.extract_and_remove_gz(gz)

    assert out == tmp_path / "file.txt"
    assert out.read_bytes() == b"hello"
    assert not gz.exists()


def test_decompress_z_rejects_bad_magic():
    with pytest.raises(ValueError):
        decompress.decompress_z(b"not a Z file")
