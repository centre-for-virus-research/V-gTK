import sys
sys.path.append("/home3/oml4h/RABV-gTK/scripts")
from CalcAlignmentCord import AlignemntCords
a = AlignemntCords()
print(a.get_products_for_range([{"start": 100, "end": 200, "product": "test_cds"}], [0, 0]))
