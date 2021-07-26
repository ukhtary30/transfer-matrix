# Tranafer Matrix method for mirror geometry

This is the repository for transfer matrix calculation of one of my paper: [Significant enhancement of light absorption in undoped graphene using dielectric multilayer system](http://aip.scitation.org/doi/abs/10.1063/1.5012604) , S. A. Nulli, **M. S. Ukhtary**, R. Saito, Appl. Phys. Lett. 112, 073101 (2018).

The summary of this paper is given by the following link of my webpages: [Cerita Riset](https://ukhtary30.github.io/significant.html)

The main unit of the program is the calculation of the absorption probability and electric field by using transfer matrix, which consists of the product between matching and propagation matrices.

The matching matrix between medium n<sub>i+1</sub> and n<sub>i</sub> is given by the following code,

    def M(x,y): #x is n_i+1, y is n_i
    A = 1 + x / y
    B = 1 - x / y
    return 0.5 * np.array([[A,B],[B,A]])

In the presence of graphene between the two media, the matching matrix is modified as,

    def MG(x,y,z): #x is n_i+1, y is n_i, z is the conductivity
    k = 2 * np.pi * y / (lmbd)
    A = 1 + x / y + k * z / (omega * eps0 * y * y)
    B = 1 - x / y + k * z / (omega * eps0 * y * y)
    C = 1 - x / y - k * z / (omega * eps0 * y * y)
    D = 1 + x / y - k * z / (omega * eps0 * y * y)
    return 0.5 * np.array([[A,B],[C,D]])

The propagation matrix in medium n<sub>i</sub> is given by the following code,

    def P(x,d): # x is the n_i, d is the thickness
    k = 2 * np.pi * x / (lmbd)
    A = 1j * k * d
    return np.array([[np.exp(-A),0],[0,np.exp(A)]])
    
    
In the paper, the thickness of each medium is the quarter of the wavelenght in the medium.



The transfer matrix of the mirror geometry with s repetition is given by the following code,

    def TF(x,s): # x is the alpha and s is the repetition
    na = x * nb
    d1 = lmbd / (4 * na)
    d2 = lmbd / (4 * nb)
    J1 = M(na,1)  @ P(na,d1) @ M(nb,na) @ P(nb,d2)
    J2 = M(na,nb) @ P(na,d1) @ M(nb,na) @ P(nb,d2)
    J3 = P(nb,d2) @ M(na,nb) @ P(na,d1) @ M(nb,na)
    J4 = P(nb,d2) @ M(na,nb) @ P(na,d1) @ M(1,na)
    JT = J1 @ matrix_power(J2,s-1) @ MG(nb,nb,sig) @ matrix_power(J3,s-1) @ J4
    return JT
    
    
