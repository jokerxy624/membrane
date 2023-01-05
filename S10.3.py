# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 22:18:17 2022

@author: joker
"""

# Simulation of dynamic evolutions for four GT-phases
# Definition and description

import numpy as np
import matplotlib.pyplot as plt

# Set of particles in random 2D Brownian motion or being attracted by large Brownian clusters if the attractive force is large enough to compete its Brownian motion

class Particle:

    def __init__(self, fixed=False, size=200, scale=0.01, ce=1e6, a=1e-6, x=None, y=None):

        self.fixed = fixed
        self.size = size
        self.scale = scale
        self.ce = ce

        if x:
            self.x = x
        else:
            self.x = np.random.uniform(0, 1)
        if y:
            self.y = y
        else:
            self.y = np.random.uniform(0, 1)

        self.cluster = None

        self.a = a

        self.moving = False
        self.target = None
        self.v = 0
        self.sin = 0
        self.cos = 0

    def step(self):

        if self.cluster:
            self.x += self.cluster.dx
            self.y += self.cluster.dy

        elif self.moving:
            x = self.x
            y = self.y
            tx = self.target.center.x
            ty = self.target.center.y

            dv = self.a * self.target.num / ((x - tx) ** 2 + (y - ty) ** 2)
            self.v += dv
            if x > tx:
                self.x -= self.v * self.sin
            else:
                self.x += self.v * self.sin
            if y > ty:
                self.y -= self.v * self.cos
            else:
                self.y += self.v * self.cos

        else:
            d = self.scale * np.random.multivariate_normal(mean=[0, 0], cov=[[1, 0], [0, 1]])

            self.x += d[0]
            self.y += d[1]

    def check_move(self, target):
        if not self.cluster:
            tx = target.center.x
            ty = target.center.y
            r = (self.x - tx) ** 2 + (self.y - ty) ** 2

            if target.num * self.size ** 2 > self.ce * r:
                self.moving = True
                self.target = target

                angle = np.arctan((tx - self.x) / (ty - self.y))
                self.sin = np.abs(np.sin(angle))
                self.cos = np.abs(np.cos(angle))

# Set of Brownian clusters in random 2D Brownian motion and will turn to Cheerios effect-driven directional attraction once the attractive force is large enough to compete with the random Brownian motion

class Cluster:

    def __init__(self, center, size=200, scale=0.01, ce=1e6, a=1e-6, stick_dis=0.02):

        self.center = center
        self.particles = [center]
        self.num = 1

        self.size = size
        self.ce = ce
        self.scale = scale
        self.a = a

        self.moving = False
        self.v = 0
        self.sin = 0
        self.cos = 0
        self.dx = 0
        self.dy = 0

        self.stick_dis = stick_dis

    def update(self, target):

        if self.moving:
            x = self.center.x
            y = self.center.y
            tx = target.center.x
            ty = target.center.y

            dv = self.a * target.num / ((x - tx) ** 2 + (y - ty) ** 2)
            self.v += dv
            if x > tx:
                self.dx = -self.v * self.sin
            else:
                self.dx = self.v * self.sin
            if y > ty:
                self.dy = -self.v * self.cos
            else:
                self.dy = self.v * self.cos

        else:
            x = self.center.x
            y = self.center.y
            tx = target.center.x
            ty = target.center.y
            r = (x - tx) ** 2 + (y - ty) ** 2

            if self.num * target.num * self.size ** 2 > self.ce * r:
                self.moving = True
                self.dx = 0
                self.dy = 0

                x = self.center.x
                y = self.center.y
                tx = target.center.x
                ty = target.center.y
                angle = np.arctan((tx - x) / (ty - y))
                self.sin = np.abs(np.sin(angle))
                self.cos = np.abs(np.cos(angle))

            else:
                d = self.scale * np.random.multivariate_normal(mean=[0, 0], cov=[[1, 0], [0, 1]]) / self.num

                self.dx = d[0]
                self.dy = d[1]

    def check_stick(self, p):

        if not p.cluster:
            x = p.x
            y = p.y

            for particle in self.particles:
                if (x - particle.x)**2 + (y - particle.y)**2 < self.stick_dis**2:
                    self.num += 1
                    self.particles.append(p)

                    p.cluster = self
                    p.fixed = True
                    return

    def check_merge(self, target):
        flag = False

        for p1 in self.particles:
            x = p1.x
            y = p1.y
            for p2 in target.particles:
                if (x - p2.x)**2 + (y - p2.y)**2 < self.stick_dis**2:
                    flag = True
                    break

        if flag:
            for p in target.particles:
                p.cluster = self

            self.particles += target.particles
            self.num += target.num

        return flag

# Initial parameters for particles
def init_particles(num_points, num_fixed, size, scale, ce, a):
    particles = [Particle(True, size, scale, ce, a, 0.5, 0.5)]
    for _ in range(1, num_fixed):
        particles.append(Particle(True, size, scale, ce, a))

    for _ in range(num_points - num_fixed):
        particles.append(Particle(False, size, scale, ce, a))

    return particles

# Parameters for Brownian clusters
def init_clusters(centers, size, scale, ce, a, stick_dis):
    clusters = []
    for center in centers:
        cluster = Cluster(center, size, scale, ce, a, stick_dis)
        clusters.append(cluster)
        center.cluster = cluster

    return clusters

# Display settings
def show(ax, particles):
    plt.pause(interval)
    ax.clear()
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.spines['top'].set_linewidth(2);
    ax.spines['bottom'].set_linewidth(2);
    ax.spines['left'].set_linewidth(2);
    ax.spines['right'].set_linewidth(2);
    plt.xticks(fontsize = 0);
    plt.yticks(fontsize = 0);
    plt.xticks([]);
    plt.yticks([]);

    x = []
    y = []
    color = []
    for particle in particles:
        x.append(particle.x)
        y.append(particle.y)
        if particle.fixed:
            color.append('blue')
        else:
            color.append('grey')

    ax.scatter(x, y, s=size, color=color, alpha=0.5)

# Rate for randomly generating new particles into the system
def add_particles(particles, num_increase, size, scale, ce, a):
    for _ in range(num_increase):
        particles.append(Particle(False, size, scale, ce, a))

# Kinematics for each particles and clusters

def step(particles, clusters):
    for i in range(1, len(particles)):
        particles[i].check_move(clusters[0])

    for i in range(1, len(clusters)):
        clusters[i].update(clusters[0])

    for particle in particles:
        particle.step()

    for cluster in clusters:
        for particle in particles:
            cluster.check_stick(particle)

    i = 1
    while i < len(clusters):
        if clusters[0].check_merge(clusters[i]):
            del clusters[i]

        else:
            i += 1

# Main program
if __name__ == '__main__':
    num_points =             # Initial numbers of isolated particles
    num_increase =           # Rate of newborn particles
    num_fixed =            # Numbers of seed particles to evolve into Brownian clusters
    stick_dis =          # Critical distance between particles to stick to each other
    scale =            # Step length of the random Brownian motion for isolated particles
    size =             # Particle size
    ce =              # Critical force between Brownian clusters/particles to activate Cheerios effect-dominated kinematics
    a =              # Factor of capillary-driven accelerated velocity
    interval =            # Interval between frames
    obs_time =            # Total frames regarded as observation duration

    particles = init_particles(num_points, num_fixed, size, scale, ce, a)
    clusters = init_clusters(particles[:num_fixed], size, scale, ce, a, stick_dis)

    fig, ax = plt.subplots(figsize=(10, 10))
    show(ax, particles)
    i = 0

    for _ in range(obs_time):

        add_particles(particles, num_increase, size, scale, ce, a)

        step(particles, clusters)

        show(ax, particles)
