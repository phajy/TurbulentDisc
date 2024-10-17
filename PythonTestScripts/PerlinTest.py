import matplotlib.pyplot as plt
from perlin_noise import PerlinNoise
from fractions import Fraction

def randomPerlin(alpha=120, beta=70, o=3):

    noise = PerlinNoise(octaves=o, seed=1)
    overlay = [[noise([i/alpha, j/beta]) for j in range(alpha)] for i in range(beta)]

    plt.imshow(overlay, cmap='gray')
    plt.xlabel('beta')
    plt.ylabel('alpha')
    plt.show()


def repeatedPerlin(alpha=120, beta=70, o=3, inverseSize = 1):

    noise = PerlinNoise(octaves=o, seed=1)

    frac = Fraction(alpha, beta)
    lim_alpha = frac.numerator * inverseSize
    lim_beta = frac.denominator * inverseSize

    overlay = [
    [
        noise([lim_alpha * i / alpha, lim_beta * j / beta], tile_sizes=[2, 3])
        for j in range(alpha)
    ]
    for i in range(beta)
]
    
    plt.imshow(overlay, cmap='gray')
    plt.xlabel('beta')
    plt.ylabel('alpha')
    plt.show()


repeatedPerlin(alpha = 240, beta = 140, inverseSize=2)
