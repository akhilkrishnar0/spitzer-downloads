from spitzer_query import SpitzerQuery


def test_field_name_instead_of_id():
    xml = b"""<?xml version="1.0"?>
<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<RESOURCE type="results">
<INFO name="QUERY_STATUS" value="OK"/>
<TABLE>
<FIELD ID="col_32" name="access_url" datatype="char" arraysize="*"/>
<DATA><TABLEDATA><TR><TD>https://example.org/file.fits</TD></TR></TABLEDATA></DATA>
</TABLE>
</RESOURCE>
</VOTABLE>"""
    table = SpitzerQuery.parse_results(xml)
    assert table.colnames == ['access_url']
    assert table['access_url'][0] == 'https://example.org/file.fits'


def test_empty_result():
    xml = b"""<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<RESOURCE type="results"><INFO name="QUERY_STATUS" value="OK"/>
<TABLE><FIELD ID="c1" name="access_url" datatype="char" arraysize="*"/>
<DATA><TABLEDATA/></DATA></TABLE></RESOURCE></VOTABLE>"""
    assert len(SpitzerQuery.parse_results(xml)) == 0
