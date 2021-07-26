# Tranafer Matrix method for mirror geometry

This is the repository for transfer matrix calculation of one of my paper: [Significant enhancement of light absorption in undoped graphene using dielectric multilayer system](http://aip.scitation.org/doi/abs/10.1063/1.5012604) , S. A. Nulli, **M. S. Ukhtary**, R. Saito, Appl. Phys. Lett. 112, 073101 (2018).

The summary of this paper is given by the following link of my webpages: [Cerita Riset](https://ukhtary30.github.io/significant.html)

The main unit of the program is the calculation of the absorption probability and electric field by using transfer matrix, which consists of the product between matching and propagation matrices.

The matching matrix between medium n<sub>i+1</sub> and n<sub>i</sub> is given by the following code,

    def M(x,y): #x is n_i+1, y is n_i
    A = 1 + x / y
    B = 1 - x / y
    return 0.5 * np.array([[A,B],[B,A]])

The propagation matrix in medium n<sub>i</sub> is given by the following code,

    def P(x,d): # x is the n_i, d is the thickness
    k = 2 * np.pi * x / (lmbd)
    A = 1j * k * d
    return np.array([[np.exp(-A),0],[0,np.exp(A)]])
    
In the paper, the thickness of each medium is the quarter of the wavelenght in the medium.


