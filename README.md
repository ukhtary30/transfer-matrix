# Tranafer Matrix method for mirror geometry

This is the repository for transfer matrix calculation of one of my paper: [Significant enhancement of light absorption in undoped graphene using dielectric multilayer system](http://aip.scitation.org/doi/abs/10.1063/1.5012604) , S. A. Nulli, **M. S. Ukhtary**, R. Saito, Appl. Phys. Lett. 112, 073101 (2018).

The summary of this paper is given by the following link of my webpages: [Cerita Riset](https://ukhtary30.github.io/significant.html)

The main unit of the program is the calculation of the absorption probability and electric field by using transfer matrix, which consists of the product between matching and propagation matrices.

The matching matrix between medium n<sub>i+1</sub> and n<sub>i</sub> is given by the following code,

    def M(x,y): #x is n_i+1, y is n_i
    A = 1 + x / y
    B = 1 - x / y
    return 0.5 * np.array([[A,B],[B,A]])

which gives the following matching matrix M<sub>i</sub>,

<img width="360" alt="Screen Shot 2021-07-27 at 9 55 36" src="https://user-images.githubusercontent.com/87349156/127078053-95a2074f-e276-4b88-8ea1-ccfe5e9ce04b.png">

In the presence of graphene between the two media, the matching matrix is modified as,

    def MG(x,y,z): #x is n_i+1, y is n_i, z is the conductivity
    k = 2 * np.pi * y / (lmbd)
    A = 1 + x / y + k * z / (omega * eps0 * y * y)
    B = 1 - x / y + k * z / (omega * eps0 * y * y)
    C = 1 - x / y - k * z / (omega * eps0 * y * y)
    D = 1 + x / y - k * z / (omega * eps0 * y * y)
    return 0.5 * np.array([[A,B],[C,D]])

which gives,

<img width="516" alt="Screen Shot 2021-07-27 at 10 00 13" src="https://user-images.githubusercontent.com/87349156/127078330-1d0e0065-b10f-49ac-bfe7-8eaae9fc93e8.png">

where we add additional term consisting of conductivity of graphene in each componenent of the matrix.

The propagation matrix in medium n<sub>i</sub> is given by the following code,

    def P(x,d): # x is the n_i, d is the thickness
    k = 2 * np.pi * x / (lmbd)
    A = 1j * k * d
    return np.array([[np.exp(-A),0],[0,np.exp(A)]])
 
 which gives the propagation matrix in medium n<sub>i</sub>, 
 
 <img width="240" alt="Screen Shot 2021-07-27 at 10 02 13" src="https://user-images.githubusercontent.com/87349156/127078478-52287170-02d4-4f5c-ad12-16f9df346871.png">
 
In the paper, the thickness of each medium (d<sub>i</sub>) is the quarter of the wavelenght in the medium.


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
    
The matrix J1 corresponds to the propagation of light when entering the mirror geometry from air to medium A and B (air --> AB). Similarly, matrix J4 corresponds to the  propagation of light leaving the geometry from medium B and to air (BA -->  air). 

Matrix J2 corresponds to the propagation of light through repeated medium AB before reaching the graphene. On the other hand, matrix J3 corresponds to the propagation of light through repeated medium BA after reaching the graphene.
