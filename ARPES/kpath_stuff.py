"""
Things I stole from berry, the module I did as TFG. It makes kpaths work pretty
seamlessly
"""
import numpy as np
import matplotlib.pyplot as plt

import itertools as it


def get_reciprocal_lattice(rLat, per=None):
    """
    Obtenemos una base de la red recíproca de nuestra red de Bravais resolviendo
    el sistema que la define, y añadiendo info extra en el caso en el que haga 
    falta.

    Parameters
    ----------
    rLat : array dimr x dimr. Contiene dimr vectores que sirvan como base de un
    espacio vectorial con dimensión dimr.

    per : list. Contiene los índices de rLat que son realmente periódicos. Si se
    deja por defecto, se considera que todos los rLat metidos son periódicos.

    Returns
    -------
    rec_lat : array dimk x dimr. Es una posible base del espacio recíproco 
    asociado al sistema
    """
    rLat = np.array(rLat)
    dimr = len(rLat[0])
    if per is None:
        dimk = dimr
        per = np.arange(dimr)
    else:
        dimk = len(per)
    
    # Generamos la matriz del sistema de ecuaciones a resolver
    if abs(np.linalg.det(rLat)) < 1e-9:
        raise Exception(f"dLat no es base de R^{dimr}! \n Hay vectores colineales en dLat")
    A = np.kron(np.eye(dimr), rLat)
    ind_b = np.arange(0, dimr**2, dimr) + np.arange(0, dimr)
    b = np.zeros(dimr**2)
    b[ind_b] = 2*np.pi
    return np.linalg.solve(A, b).reshape(dimr, dimr)[per]

def brillouin_zone_structure(recLat, per=None):
    if per is None:
        dimk = recLat.shape[0]
        per = np.arange(0, dimk, 1)
    else:
        dimk = len(per)

    # We create the mesh and eliminate the value for zero
    Nrange = np.arange(-1, 2, 1)
    Nmesh = np.array(np.meshgrid(*(dimk*[Nrange])))[::-1].reshape(dimk, 3**dimk)
    Nmesh = np.delete(Nmesh, 3**dimk//2, axis=1)
    
    Gmesh = np.einsum("i...,ir->...r", Nmesh, recLat[:,per])
    Gmod = np.linalg.norm(Gmesh, axis=1)

    per_list = list(it.combinations(tuple(range(len(Gmod))), dimk))
    sols = []
    ind_sols = []
    for perm in per_list:
        try:
            sols.append(np.linalg.solve(Gmesh[perm,:], Gmod[list(perm)]**2/2))

        except np.linalg.LinAlgError:
            pass
    sols = np.array(sols)
    mod_sols = np.sqrt(np.sum(sols**2, axis=1))
    # TODO: relate the limit to the dimension of the space and the module of the 
    # lattice parameter
    dots = sols @ Gmesh.T                      # (nsol, nG)
    rhs = (Gmod**2)/2                          # (nG,)
    mask = np.all(dots <= rhs + 1e-10, axis=1) # boolean mask
    sols_f1 = sols[mask]
    #hull = ConvexHull(sols_f1)
    #vertices = hull.vertices 
    return sols_f1


def rec_lattice(rLat, per=None):
    recLat = get_reciprocal_lattice(rLat, per)
    BZ_points = brillouin_zone_structure(recLat, per)
    return recLat, BZ_points


def path_chart(points, nk, recLat, mode=0):
    """
    Creamos un camino en un espacio de parámetros con una determinada densidad
    de puntos en cada camino. Por comodidad, los ptos están entre 0 y 1
    
    Los puntos NO eran equiespaciados en el espacio recíproco!!
    
    Hay 2 opciones para el input: 
        -fijamos el n del primer segmento y con él, fijamos dk
        -fijamos el n de cada segmento. Los dk son diferentes!!
    
    Parameters
    ----------
    points : Lista. Contiene puntos en dimk

    Returns
    -------
    kpath : array dimk x nk?. Contiene todos los ptos del espacio recíproco por 
    los que pasa el camino

    kind : array dimk x nk?. Contiene los índices en los que hay un cambio de 
    segmento para poderlos graficar de forma eficiente

    kdist : array dimk x nk?. Contiene la distancia recorrida en el esp. rec.
    para graficar las bandas acordemente a la misma
    """
    # Pasamos los puntos a array para trabajar mejor con ellos
    points = np.array(points)

    # Vectores que unen los ptos en el espacio recíproco y su módulo
    dPoints = points[1:] - points[:-1]

    # Calculamos nuestros ptos en la red recíproca
    k_mod = np.sum(np.abs(dPoints)**2, axis=1)
    # -------------------------------------------------------------------------
    # SI nk ES UN ENTERO
    if mode == 0:
        # Calculamos el dk global, e intentamos hacer que sea parecido a este 
        # en el resto de puntos
        dk = k_mod[0]/nk
        # Guardamos los segmentos en una lista que después concatenamos
        l_seg = [np.linspace(points[0], points[1], nk)]
        # Guardamos los puntos cuando comienza un segmento nuevo
        kind = [0, nk]
        # Añadimos los que quedan
        for i in range(1, len(points)-1):
            naux = int(k_mod[i]/dk)
            l_seg.append(np.linspace(points[i], points[i+1], naux))
            kind.append(naux)
        # Concatenamos todos los resultados en uno solo
        kpath = np.concatenate(tuple(l_seg), axis = 0)
    # -------------------------------------------------------------------------
    # SI nk ES UNA LISTA O ARRAY
    if mode == 1:
        # Global dk for all segments
        dk_g = np.sum(k_mod)/nk

        # First segment subdivision
        nk_parc = int(k_mod[0]/dk_g)
        l_seg = [np.linspace(points[0], points[1], nk_parc)]
        kiond = [0, nk_parc]
        for i in range(1, len(points)-1):
            nk_parc = int(k_mod[i]/dk_g)
            l_seg.append(np.linspace(points[i], points[i+1], naux))
            kind.append(nk_parc)
        kpath = np.concatenate(tuple(l_seg), axis = 0)

    if mode == 2:
        l_seg = []
        kind = [0]
        # Hacemos y encadenamos los linspace
        for i in range(len(points) - 1):
            l_seg.append(np.linspace(points[i], points[i+1], nk[i]))
            kind.append(nk[i])
        kpath = np.concatenate(tuple(l_seg), axis = 0)
    # -------------------------------------------------------------------------
    # Calculamos las diferentes distancias de nuestros ks
    kdist = np.sqrt(np.sum((kpath[1:] - kpath[:-1])**2, axis=1))
    
    kind[-1] = kind[-1] -1
    # Devolvemos las distancias acumuladas para poder graficar con sentido
    return kpath.T, np.cumsum(kind), np.concatenate((np.array([0]), np.cumsum(kdist)), axis=0)


def plot_1BZ(recLat, BZ_points):
    # Points of the reciprocal lattice
    # TODO: Make it so it always ends the square nicely with convex or something
    relevant_pts = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])
    Lat_points = relevant_pts @ recLat
    fig, ax = plt.subplots(dpi=200)
    ax.plot(Lat_points[:,0], Lat_points[:,1], c='gray', marker='.')
    ax.scatter(BZ_points[:,0], BZ_points[:,1], color='blue')
    ax.set_xlabel('$k_x$')
    ax.set_ylabel('$k_y$')
    return fig, ax