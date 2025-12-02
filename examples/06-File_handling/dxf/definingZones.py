from hughes2d import NonConvexDomain

filename = "examples/06-File_handling/dxf/config_simple_zones.dxf"

MyDomain = NonConvexDomain()
MyDomain.importFromDXF(filename)
MyDomain.show()
