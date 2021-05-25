import math_tools as mt
import classes as c
import tools as ts
import sel as s

def app():

    localKs = []
    localbs = []
    Kv = []
    b = []
    Tv = []

    print("IMPLEMENTACION DEL METODO DE LOS ELEMENTOS FINITOS\n")
    print("\t- TRANSFERENCIA DE CALOR\n" )
    print("\t- 2 DIMENSIONES\n")
    print("\t- FUNCIONES DE FORMA LINEALES\n")
    print("\t- PESOS DE GALERKIN\n")
    print("\t- MALLA TRIANGULAR IRREGULAR\n")
    print("*********************************************************************************\n\n")
    
    m = c.Mesh()
    filename = "test2"

    ts.leerMallayCondiciones(m, filename)
    print("Datos obtenidos correctamente\n********************\n")

    s.crearSistemasLocales(m, localKs, localbs)
    print("******************************\n")

    mt.ZeroesFirst(Kv, m.getSize(c.Sizes.NODES.value))
    mt.ZeroesThird(b, m.getSize(c.Sizes.NODES.value))

    s.ensamblaje(m, localKs, localbs, Kv, b)
    print("******************************\n")

    s.applyNeumann(m, b)
    print("******************************\n")

    s.applyDirichlet(m, Kv, b)
    print("******************************\n")

    mt.ZeroesThird(Tv, len(b))
    s.calculate(Kv, b, Tv)

    ts.writeResults(m, Tv, filename)

app()