cimport cython

@cython.cdivision(True)
cpdef double saturation_pressure_water(double temperature):
    # IAPWS97 - Region 4
    cdef double n1 = 0.11670521452767E4
    cdef double n2 = -0.72421316703206E6
    cdef double n3 = -0.17073846940092E2
    cdef double n4 = 0.12020824702470E5
    cdef double n5 = -0.32325550322333E7
    cdef double n6 = 0.14915108613530E2
    cdef double n7 = -0.48232657361591E4
    cdef double n8 = 0.40511340542057E6
    cdef double n9 = -0.23855557567849
    cdef double n10 = 0.65017534844798E3
    cdef double pstar = 1e6
    cdef double Tstar = 1
    cdef double theta, A, B, C, ps

    theta = (temperature / Tstar) + (n9 / ((temperature / Tstar) - n10))
    A = theta ** 2 + n1 * theta + n2
    B = n3 * theta ** 2 + n4 * theta + n5
    C = n6 * theta ** 2 + n7 * theta + n8

    ps = pstar * ((2 * C) / (-B + (B ** 2 - 4 * A * C) ** 0.5)) ** 4
    return ps

@cython.cdivision(True)
cpdef double saturation_Temperature_water(double pressure):
    # IAPWS97 - Region 4
    cdef double n1 = 0.11670521452767E4
    cdef double n2 = -0.72421316703206E6
    cdef double n3 = -0.17073846940092E2
    cdef double n4 = 0.12020824702470E5
    cdef double n5 = -0.32325550322333E7
    cdef double n6 = 0.14915108613530E2
    cdef double n7 = -0.48232657361591E4
    cdef double n8 = 0.40511340542057E6
    cdef double n9 = -0.23855557567849
    cdef double n10 = 0.65017534844798E3
    cdef double pstar = 1e6
    cdef double Tstar = 1
    cdef double beta, D, E, F, G, Ts

    beta = (pressure / pstar) ** 0.25
    E = beta ** 2 + n3 * beta + n6
    F = n1 * beta ** 2 + n4 * beta + n7
    G = n2 * beta ** 2 + n5 * beta + n8
    D = (2 * G) / (-F - (F ** 2 - 4 * E * G) ** 0.5)

    Ts = Tstar * 0.5 * (n10 + D - ((n10 + D) ** 2 - 4 * (n9 + n10 * D)) ** 0.5)
    return Ts