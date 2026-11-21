"""Sistema acoplado fluido–sólido elástico en 3D (Fase 6b).

Es el paso que falta para poder comparar el FEM contra **Faran** (D45): hasta acá el 3D sólo sabía de
fluidos (Helmholtz escalar) y de un dispersor con contraste de impedancia, no de un sólido que además
sostiene ondas de corte.

Incógnitas
----------
- ``u`` : desplazamiento del sólido, 3 grados de libertad por nodo (``elastic3d``).
- ``p_s`` : presión **dispersada** en el fluido, 1 por nodo.

Forma débil
-----------
En el sólido, con la tracción del fluido como condición de borde:

    ∫σ:ε(v) − ω²∫ρ u·v + ∫_Γ p_s (n_s·v) dS = −∫_Γ p_i (n_s·v) dS

En el fluido (cáscara ``a < r < b``), con ABC en ``r = b``:

    ∫(1/ρ)∇p_s·∇q − (k²/ρ)∫p_s q − (ik/ρ)∫_∂Ω p_s q + ω²∫_Γ q (u·n_s) dS
        = (1/ρ)∫_Γ q ∂p_i/∂n_s dS

Los dos signos que hay que cuidar y que están verificados en los tests:

- ``n_s`` es la normal **saliente del sólido** (``interface_faces`` se orienta así al construir la
  malla), y en el fluido la normal saliente del dominio es ``−n_s``; de ahí el signo del término ABC
  del lado del fluido y el ``+ω²`` del acoplamiento.
- La continuidad es de **desplazamiento** (``u_f·n = u_s·n``) y de **tracción** (``σ·n_s = −p n_s``),
  no de presión: el sólido no tiene presión.

Sistema resultante (no simétrico tal cual, y está bien):

    [ K − ω²M      C  ] [u ]   [ −∫p_i (n_s·v)          ]
    [ ω²Cᵀ        A_ff] [p_s] = [ (1/ρ)∫q ∂p_i/∂n_s      ]

``C`` es la matriz de masa de la superficie de interfaz proyectada sobre la normal. No se usa
multiplicador de Lagrange ni proyección: la malla comparte nodos en la interfaz.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..materials import ElasticSolid, Material
from ..meshing3d_sphere import SphereShellMesh
from . import elastic3d, helmholtz3d

__all__ = ["matrices_acopladas", "resolver", "comparar_con_faran"]


def _matriz_masa_cara(nodes: np.ndarray, faces: np.ndarray) -> sp.csr_matrix:
    """Matriz de masa ``∫NᵢNⱼ`` sobre triángulos P1, ensamblada sobre los nodos de ``faces``."""
    p = nodes[faces]
    areas = 0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1)
    local = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]) / 12.0
    filas = np.repeat(faces, 3, axis=1).reshape(-1)
    cols = np.tile(faces, (1, 3)).reshape(-1)
    vals = (areas[:, None, None] * local[None, :, :]).reshape(-1)
    n = nodes.shape[0]
    return sp.coo_matrix((vals, (filas, cols)), shape=(n, n)).tocsr()


def _matriz_fluido(
    mesh: SphereShellMesh, fluid: Material, frequency: float, abc_order: int = 1
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``A_ff`` del fluido en su **espacio reducido**, y las caras ABC con sus normales.

    El espacio reducido importa: si se arma sobre todos los nodos, los que viven sólo dentro de la
    esfera quedan con fila vacía (medido: 44 filas) y el sistema global sale **singular**. El error
    que da SuperLU es "Factor is exactly singular", que no dice nada sobre la causa.
    """
    nodos_fluido = np.unique(mesh.tets_fluid)
    mapa = np.full(mesh.n_nodes, -1, dtype=np.int64)
    mapa[nodos_fluido] = np.arange(nodos_fluido.size)
    nodes = mesh.nodes[nodos_fluido]
    tets = mapa[mesh.tets_fluid]
    k = complex(fluid.wavenumber(frequency))
    rho = float(fluid.rho)

    grads, vols = helmholtz3d.tet_gradients(nodes, tets)
    # rigidez: ∫(1/ρ)∇N·∇N ; masa: ∫(k²/ρ)NN  (P1 tet: ∫NN = V/20·(2 diag, 1 fuera))
    # suma sobre el índice ESPACIAL (d), no sobre los nodos: se quiere la matriz 4×4 nodo-nodo
    ke = np.einsum("mid,mjd->mij", grads, grads) * vols[:, None, None] / rho
    masa_local = np.full((4, 4), 1.0)
    np.fill_diagonal(masa_local, 2.0)
    me = (k**2 / rho) * vols[:, None, None] * masa_local[None, :, :] / 20.0
    n = nodes.shape[0]
    dofs = tets
    filas = np.repeat(dofs, 4, axis=1).reshape(-1)
    cols = np.tile(dofs, (1, 4)).reshape(-1)
    a = sp.coo_matrix(((ke - me).reshape(-1), (filas, cols)), shape=(n, n)).tocsr()

    # caras de frontera del fluido: la esfera exterior (ABC) y la interfaz (que va aparte)
    caras = helmholtz3d_faces(tets)
    clave_interfaz = {tuple(sorted(f)) for f in mesh.interface_faces}
    exterior = np.array([f for f in caras if tuple(sorted(f)) not in clave_interfaz], dtype=np.int64)
    if exterior.size == 0:
        raise RuntimeError("el fluido no tiene frontera exterior: revisar la malla")

    p = nodes[exterior]
    normales = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    normales /= np.linalg.norm(normales, axis=1)[:, None]
    # la normal de boundary_faces ya es saliente del fluido; se verifica con el centroide
    afuera = np.einsum("ij,ij->i", normales, p.mean(axis=1)) > 0
    if not np.all(afuera):
        raise RuntimeError("hay caras exteriores del fluido con normal entrante")

    # ABC sobre la esfera exterior. El orden 2 (BGT-2) no es un lujo: con k₀b ≈ 2.5 el ABC de
    # impedancia plana (orden 1) deja pasar buena parte de la onda, y el error que se mide contra
    # Faran es del ABC, no de la discretización (medido: 42 % con orden 1).
    mid = nodes[exterior].mean(axis=1)
    r_cara = np.linalg.norm(mid, axis=1)
    c_f, d_f = helmholtz3d.abc_coefficients(k, r_cara, order=abc_order)
    # con orden 1 los coeficientes salen escalares (no usan el radio): se replican por cara
    c_f = np.broadcast_to(np.asarray(c_f, dtype=complex), r_cara.shape)
    d_f = np.broadcast_to(np.asarray(d_f, dtype=complex), r_cara.shape)
    masa_local = helmholtz3d.face_mass(nodes, exterior)          # (F, 3, 3)
    w = 1.0 / rho
    a = a + _acumular(nodes.shape[0], exterior, (-c_f * w)[:, None, None] * masa_local)
    if abc_order == 2:
        rig_local = helmholtz3d.surface_stiffness(nodes, exterior)
        a = a + _acumular(nodes.shape[0], exterior, (d_f * w)[:, None, None] * rig_local)
    return a.tocsr(), nodos_fluido, mapa, exterior, normales


def _acumular(n: int, faces: np.ndarray, locales: np.ndarray) -> sp.csr_matrix:
    """Acumula matrices locales ``(F, 3, 3)`` de las caras en una matriz global ``(n, n)``."""
    filas = np.repeat(faces, 3, axis=1).reshape(-1)
    cols = np.tile(faces, (1, 3)).reshape(-1)
    return sp.coo_matrix((np.asarray(locales).reshape(-1), (filas, cols)), shape=(n, n)).tocsr()


def helmholtz3d_faces(tets: np.ndarray) -> np.ndarray:
    """Caras de frontera de una región, reutilizando la utilidad de `meshing3d`."""
    from ..meshing3d import boundary_faces

    return boundary_faces(tets)


def matrices_acopladas(
    mesh: SphereShellMesh,
    fluid: Material,
    solid: ElasticSolid,
    frequency: float,
    abc_order: int = 1,
) -> dict[str, Any]:
    """Ensambla los cuatro bloques y devuelve todo lo necesario para resolver."""
    a_ff, nodos_fluido, mapa_f, exterior, normales_ext = _matriz_fluido(
        mesh, fluid, frequency, abc_order
    )

    # El sólido vive en un espacio REDUCIDO: sólo los nodos que tocan sus tetraedros. Si se arma la
    # rigidez sobre todos los nodos, las filas de los nodos sólo-fluido quedan en cero y el sistema
    # global sale singular (y el error que da es de dimensiones, no de física).
    nodos_solidos = np.unique(mesh.tets_solid)
    mapa = np.full(mesh.n_nodes, -1, dtype=np.int64)
    mapa[nodos_solidos] = np.arange(nodos_solidos.size)
    n_solid = nodos_solidos.size
    k, m = elastic3d.p1_elastic3d_matrices(
        mesh.nodes[nodos_solidos], mapa[mesh.tets_solid], solid
    )

    # acoplamiento: C[nodo_fluido, dof_solido] = ∫ N_i N_j (n_s · e_dir) dS
    caras = mesh.interface_faces
    normales = mesh.interface_normals()
    c = sp.lil_matrix((nodos_fluido.size, 3 * n_solid), dtype=complex)
    p = mesh.nodes[caras]
    areas = mesh.interface_areas()
    local = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]) / 12.0
    for f in range(caras.shape[0]):
        for i in range(3):
            for j in range(3):
                peso = areas[f] * local[i, j]
                nodo_s = mapa[caras[f, j]]
                if nodo_s < 0:
                    continue
                fila = mapa_f[caras[f, i]]
                if fila < 0:
                    continue
                for d in range(3):
                    c[fila, 3 * nodo_s + d] += peso * normales[f, d]
    c = c.tocsr()

    # el bloque de acoplamiento del sólido usa los mismos pesos: Cᵀ, y el del fluido ω²·Cᵀ
    return {
        "K": k, "M": m, "A_ff": a_ff, "C": c, "mapa": mapa, "mapa_f": mapa_f,
        "nodos_fluido": nodos_fluido, "nodos_solidos": nodos_solidos,
        "n_solid": n_solid, "exterior": exterior, "normales_exterior": normales_ext,
        "omega": 2.0 * np.pi * frequency, "k0": complex(fluid.wavenumber(frequency)),
    }


def resolver(
    mesh: SphereShellMesh,
    fluid: Material,
    solid: ElasticSolid,
    frequency: float,
    onda: Any = None,
    abc_order: int = 1,
    fijar_rigidos: bool | None = None,
    formulacion: str = "total",
) -> dict[str, Any]:
    """Resuelve la esfera elástica sumergida. Dos formulaciones, y la elección no es cosmética.

    ``formulacion="total"`` (por defecto) resuelve el desplazamiento completo y, con
    ``fijar_rigidos=True``, impone traslación media nula para domar los seis modos rígidos. Es la que
    reproduce Faran al **4.29 %** en metales. Pero la condición **suprime una respuesta dipolar
    física** (Feuillade y Clay 1999, §IV; D48), así que ese 4 % es en parte una cancelación
    afortunada, y con materiales livianos el error crece como ``3/(2g+1)``: 13 % en hueso, catastrófico
    en tejido.

    ``formulacion="dispersada"`` resuelve ``u_s = u − u_i`` con el **operador de contraste** aplicado
    a la incidente. Es la correcta —el caso transparente da ``|p_s| = 1e-4`` sin restringir nada— pero
    numéricamente es una **resta de términos grandes**: para un metal ``u_s ≈ −u_i`` con la rigidez del
    acero, y cualquier inconsistencia de la incidente se amplifica. Medido contra Faran: 19.49 % en
    acero y 44.14 % en aluminio. No se usa como default hasta entender esa amplificación (D49).
    """
    """Resuelve el sistema acoplado para una onda plana incidente.

    Parámetros
    ----------
    onda : `helmholtz3d.PlaneWave3D` o ``None`` (por defecto, propagación en ``+z``).
    """
    if onda is None:
        onda = helmholtz3d.plane_wave(frequency, fluid, np.array([0.0, 0.0, 1.0]))

    if formulacion not in ("dispersada", "total"):
        raise ValueError(f"formulacion tiene que ser 'dispersada' o 'total'; llegó {formulacion!r}")
    if fijar_rigidos is None:
        # automático: la restricción sólo tiene sentido en la formulación total. Dejarlo en un
        # default booleano hacía que pedir la dispersada reventara sin decir por qué.
        fijar_rigidos = formulacion == "total"

    datos = matrices_acopladas(mesh, fluid, solid, frequency, abc_order)
    omega = datos["omega"]
    k_solid, m_solid, a_ff, c = datos["K"], datos["M"], datos["A_ff"], datos["C"]
    mapa, nodos_solidos, n_solid = datos["mapa"], datos["nodos_solidos"], datos["n_solid"]
    nodos_fluido = datos["nodos_fluido"]

    caras = mesh.interface_faces
    normales = mesh.interface_normals()
    p_cara = mesh.nodes[caras]
    areas = mesh.interface_areas()
    masa_peso = np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]]) / 12.0

    # término independiente del sólido: −∫ p_i (n_s·v)
    b_s = np.zeros(3 * n_solid, dtype=complex)
    # término independiente del fluido: (1/ρ)∫ q (∂p_i/∂n_s)
    mapa_f = datos["mapa_f"]
    b_f = np.zeros(datos["nodos_fluido"].size, dtype=complex)
    rho = float(fluid.rho)
    for f in range(caras.shape[0]):
        puntos = p_cara[f]
        p_i = onda.value(puntos)
        dp_i = onda.normal_derivative(puntos, np.tile(normales[f], (3, 1)))
        for i in range(3):
            nodo_f = mapa_f[caras[f, i]]
            for j in range(3):
                peso = areas[f] * masa_peso[i, j]
                nodo_s = mapa[caras[f, j]]
                if nodo_s >= 0:
                    b_s[3 * nodo_s:3 * nodo_s + 3] -= peso * p_i[j] * normales[f]
                b_f[nodo_f] += peso / rho * dp_i[j]

    # bloques globales: [u (3N) ; p_s (N)]
    n = mesh.n_nodes
    cero = sp.csr_matrix((3 * n_solid, 3 * n_solid))
    n_total = 3 * n_solid + nodos_fluido.size
    a_global = sp.bmat([
        [k_solid - omega**2 * m_solid, c.T],
        [omega**2 * c, a_ff],
    ], format="csc")

    # Desplazamiento de la onda incidente prolongado al interior del sólido. Se usa la relación de
    # momento del fluido, `u = ∇p/(ρ₀ω²)`, que es la misma con la que se impone la continuidad: eso es
    # lo que hace que el caso transparente cierre sin restricción.
    u_inc = np.asarray(onda.gradient(mesh.nodes[nodos_solidos]), dtype=complex) / (
        float(fluid.rho) * omega**2
    )
    u_inc_vec = u_inc.reshape(-1)

    if formulacion == "dispersada":
        # Sólido en campo DISPERSADO: u = u_i + u_s. La ecuación del sólido queda
        #   (K − ω²M)·u_s = −[(K − ω²M)·u_i + ∫p_i(n·v)] − Cᵀ·p_s
        # y el corchete es el residuo del operador ELÁSTICO aplicado a la incidente. Para un cuerpo
        # indistinguible del fluido ese residuo se anula por consistencia (el término de volumen de
        # ∫σ:ε(v) cancela contra el de borde), así que u_s = 0 sin imponer nada.
        # El fluido pierde su dato de Neumann: ρ₀ω²(u_i·n) y ∂p_i/∂n son el mismo número, y se
        # cancelan. Todo lo que entra al sistema entra por el residuo del sólido.
        # Término de contraste, no residuo del operador completo. La identidad
        #   (K_fluido − ω²M_fluido)·u_i = −∮p_i(n·v)
        # se cumple en el continuo, pero al restarla numéricamente queda un residuo O(h²) (medido:
        # 0.35 % de la tracción) que **empuja un campo dispersado espurio** de ~0.06 Pa, del orden del
        # 18 % del campo real de un metal.
        #
        # La forma correcta —la misma lección que D41 en el caso escalar— es aplicar sólo el CONTRASTE
        # del operador sobre la incidente. Con el operador del fluido armado sobre la MISMA malla y la
        # MISMA cuadratura, un cuerpo indistinguible del fluido da matrices idénticas y el término es
        # **exactamente cero**, sin depender de h.
        lam_f = float(fluid.rho) * float(fluid.c) ** 2
        k_fluido, m_fluido = elastic3d.p1_elastic3d_matrices_constantes(
            mesh.nodes[nodos_solidos], mapa[mesh.tets_solid], lam_f, 0.0, float(fluid.rho)
        )
        contraste = (k_solid - omega**2 * m_solid) - (k_fluido - omega**2 * m_fluido)
        residuo_incidente = contraste @ u_inc_vec
        b_global = np.concatenate([-residuo_incidente, np.zeros(nodos_fluido.size, dtype=complex)])
    else:
        b_global = np.concatenate([b_s, b_f])

    multiplicadores = None

    if fijar_rigidos and formulacion == "dispersada":
        raise ValueError(
            "fijar_rigidos es de la formulación 'total'. En 'dispersada' no hace falta: el "
            "desplazamiento medio de u_s ya es el correcto (y restringirlo suprimiría el dipolo)."
        )

    if fijar_rigidos:
        # El sólido libre tiene seis modos de cuerpo rígido SIN rigidez elástica: la única que los
        # contiene es la reacción del fluido, que es chica. Con una onda plana la fuerza neta exacta
        # es CERO (por simetría del problema), pero la malla no la cancela exactamente, y ese residuo
        # excita la traslación sin freno: medido, el desplazamiento normal en la interfaz salía 20×
        # el exacto y la amplitud dispersada 16 % baja.
        #
        # Se impone la traslación media nula como restricción de Lagrange (válido para una onda plana
        # sobre un cuerpo simétrico) y se DEVUELVEN los multiplicadores: son la fuerza y el torque
        # netos que la malla no supo cancelar. Así el desbalance queda medido en vez de escondido
        # dentro de un desplazamiento espurio.
        modos = elastic3d.rigid_body_modes(mesh.nodes[nodos_solidos])
        g = sp.csr_matrix(modos.T)                       # (6, 3N)
        cero6 = sp.csr_matrix((6, 6))
        a_global = sp.bmat([
            [a_global, sp.vstack([g.T, sp.csr_matrix((nodos_fluido.size, 6))])],
            [sp.hstack([g, sp.csr_matrix((6, nodos_fluido.size))]), cero6],
        ], format="csc")
        b_global = np.concatenate([b_global, np.zeros(6, dtype=complex)])

    lu = spla.splu(a_global)
    sol = lu.solve(b_global)
    x = sol[:n_total]
    if fijar_rigidos:
        multiplicadores = sol[n_total:]
    u = x[: 3 * n_solid].reshape(n_solid, 3)
    p_s_reducido = x[3 * n_solid:]

    # se devuelve en el espacio COMPLETO de nodos para que la interpolación y las gráficas no
    # tengan que saber de las reducciones internas. En la formulación dispersada lo resuelto es
    # `u_s`, así que el desplazamiento TOTAL del sólido es `u_i + u_s`.
    desplazamiento = np.zeros((n, 3), dtype=complex)
    desplazamiento[nodos_solidos] = u + (u_inc if formulacion == "dispersada" else 0.0)
    desplazamiento_dispersado = np.zeros((n, 3), dtype=complex)
    desplazamiento_dispersado[nodos_solidos] = u
    p_s = np.zeros(n, dtype=complex)
    p_s[datos["nodos_fluido"]] = p_s_reducido
    return {
        "u": desplazamiento, "p_s": p_s, "presion_total": p_s + onda.value(mesh.nodes),
        "frequency": frequency, "onda": onda, "datos": datos, "n_solid": n_solid,
        "u_dispersado": desplazamiento_dispersado, "u_incidente": u_inc,
        "formulacion": formulacion,
        # fuerza y torque netos que la malla no canceló (None si no se impuso la restricción)
        "residuo_rigido": multiplicadores,
    }


def comparar_con_faran(
    resultado: dict[str, Any], mesh: SphereShellMesh, fluid: Material, solid: ElasticSolid,
    puntos: np.ndarray, direccion=(0.0, 0.0, 1.0),
) -> dict[str, Any]:
    """Contrasta la presión dispersada del FEM contra la serie de Faran en los mismos puntos."""
    from ..analytic import faran

    coef = faran.coeficientes(resultado["frequency"], fluid, solid, mesh.radius)
    exacto = faran.campo_dispersado(coef, puntos, mesh.radius, direccion)

    # interpolación P1 del campo del FEM en los puntos pedidos: se usa el tet que contiene el punto
    fem = _interpolar(mesh, resultado["p_s"], puntos)
    peso = np.abs(exacto)
    util = peso > 1e-3 * peso.max()
    error = np.abs(fem[util] - exacto[util]).max() / peso[util].max() if np.any(util) else np.inf
    return {
        "puntos": puntos, "fem": fem, "exacto": exacto, "error_relativo": float(error),
        "amplitud_maxima": float(np.abs(exacto).max()), "ka": coef["ka"],
        "n_puntos": int(util.sum()),
    }


def _interpolar(
    mesh: SphereShellMesh, valores: np.ndarray, puntos: np.ndarray, vecinos: int = 32
) -> np.ndarray:
    """Interpolación P1 de un campo nodal en puntos arbitrarios, buscando el tet contenedor.

    Si un punto cae fuera de la malla se devuelve ``nan`` y el llamador lo descarta en vez de inventar
    un valor (los puntos dentro de la esfera salen ``nan``, que es lo correcto: ahí no hay fluido).

    Cada punto se busca sólo entre los ``vecinos`` tetraedros con centroide más cercano, con un
    KD-tree. Probarlos todos sería correcto pero **8× más lento**: con una grilla de 81×81 para el
    mapa, 6 561 puntos × 8 672 tets pasaba de 15 s, y el solve entero son 1.8 s. El radio de
    búsqueda está de sobra porque los centroides de una malla de tamaño ``h`` están a menos de ``h``
    entre sí.
    """
    from scipy.spatial import cKDTree

    puntos = np.atleast_2d(np.asarray(puntos, dtype=float))
    salida = np.full(puntos.shape[0], np.nan + 0j)
    tets = mesh.tets_fluid
    p = mesh.nodes[tets]                       # (M,4,3)
    v0 = p[:, 1] - p[:, 0]
    v1 = p[:, 2] - p[:, 0]
    v2 = p[:, 3] - p[:, 0]
    det = np.einsum("ij,ij->i", v0, np.cross(v1, v2))
    if det.size == 0:
        return salida

    arbol = cKDTree(p.mean(axis=1))
    k = int(min(max(vecinos, 4), det.size))
    _, candidatos = arbol.query(puntos, k=k)
    candidatos = np.atleast_2d(candidatos)

    for i, x in enumerate(puntos):
        cand = candidatos[i]
        d = x - p[cand, 0]
        l1 = np.einsum("ij,ij->i", d, np.cross(v1[cand], v2[cand])) / det[cand]
        l2 = np.einsum("ij,ij->i", v0[cand], np.cross(d, v2[cand])) / det[cand]
        l3 = np.einsum("ij,ij->i", v0[cand], np.cross(v1[cand], d)) / det[cand]
        l0 = 1.0 - l1 - l2 - l3
        dentro = (l0 >= -1e-9) & (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
        if not np.any(dentro):
            continue
        e = int(cand[np.flatnonzero(dentro)[0]])
        j = int(np.flatnonzero(dentro)[0])
        salida[i] = (l0[j] * valores[tets[e, 0]] + l1[j] * valores[tets[e, 1]]
                     + l2[j] * valores[tets[e, 2]] + l3[j] * valores[tets[e, 3]])
    return salida
