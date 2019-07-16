#!/usr/bin/python3

import numpy as np

NEWTON = True


def acceleration(star, bh1, bh2, **kwargs):
    if 'G' in kwargs:
        G = kwargs['G']
    else:
        G = 1.0

    if 'c_light' in kwargs:
        c_light = kwargs['c_light']
    else:
        c_light = 1.0

    if 'mass' in kwargs:
        mass = kwargs['mass']
    else:
        mass = 1.0

    if 'massratio' in kwargs:
        q = kwargs['massratio']
    else:
        q = 1.0

    dist_star_bh1 = star.pos - bh1.pos
    dist_star_bh2 = star.pos - bh2.pos

    acc_from_bh1 = (dist_star_bh1) * (-1 * bh1.mass * G) / \
        (np.linalg.norm(dist_star_bh1) ** 3)
    acc_from_bh2 = (dist_star_bh2) * (-1 * bh2.mass * G) / \
        (np.linalg.norm(dist_star_bh2) ** 3)

    if NEWTON:
        return acc_from_bh1 + acc_from_bh2

    second_acc_from_bh1 = c_light ** -2 * (G * bh1.mass * dist_star_bh1) / (np.linalg.norm(dist_star_bh1) ** 2) * \
        np.linalg.norm(star.vel) ** 2 + 2 * np.linalg.norm(bh1.vel) ** 2 - 4 * np.dot(star.vel, bh1.vel) - 3/2 * np.dot(dist_star_bh1, bh1.vel) ** 2 \
        - 4 * ((G * bh2.mass / dist_star_bh2) + (G * bh1.mass / dist_star_bh1)) - \
        (G * bh2.mass / dist_star_bh2) + 0.5 * np.dot(dist_star_bh1, bh1.acc)

    second_acc_from_bh2 = c_light ** -2 * (G * bh2.mass * dist_star_bh2) / (np.linalg.norm(dist_star_bh2) ** 2) * \
        np.linalg.norm(star.vel) ** 2 + 2 * np.linalg.norm(bh2.vel) ** 2 - 4 * np.dot(star.vel, bh2.vel) - 3/2 * np.dot(dist_star_bh2, bh2.vel) ** 2 \
        - 4 * ((G * bh1.mass / dist_star_bh1) + (G * bh2.mass / dist_star_bh2)) - \
        (G * bh1.mass / dist_star_bh1) + 0.5 * np.dot(dist_star_bh2, bh2.acc)

    # return second_acc_from_bh1 + second_acc_from_bh2

    third_acc_from_bh1 = c_light ** -2 * G * bh1.mass / \
        (np.linalg.norm(dist_star_bh1) ** 2) * np.dot(dist_star_bh1,
                                                      4 * star.vel - 3 * bh1.vel) * (star.vel - bh1.vel)

    third_acc_from_bh2 = c_light ** -2 * G * bh2.mass / \
        (np.linalg.norm(dist_star_bh2) ** 2) * np.dot(dist_star_bh2,
                                                      4 * star.vel - 3 * bh2.vel) * (star.vel - bh2.vel)

    fourth_acc_from_bh1 = 7 / (2 * c_light ** 2) * G * bh1.mass * bh1.acc / \
        np.linalg.norm(dist_star_bh1) + 0.0  # less significant terms

    fourth_acc_from_bh2 = 7 / (2 * c_light ** 2) * G * bh2.mass * bh2.acc / \
        np.linalg.norm(dist_star_bh1) + 0.0  # less significant terms

    # done with terms, now sum it all

    acc_from_bh1 += second_acc_from_bh1 + third_acc_from_bh1 + fourth_acc_from_bh1
    acc_from_bh2 += second_acc_from_bh2 + third_acc_from_bh2 + fourth_acc_from_bh2

    return acc_from_bh1 + acc_from_bh2


if __name__ == "__main__":
    pass
