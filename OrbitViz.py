from vpython import * #Already imports numpy and presumably, math

def solvekep(M): #Iterative method to determine the solution to Kepler's equation
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M #Should this line be 
    while abs(Mguess - M) > 0.0001:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (Mtrue-(Eguess-e*sin(Eguess)))/ (e*cos(Eguess)-1)
    return Eguess

def locateAsteroid(wprime, Oprime, iprime, a, e, E, r):
    r.x = (a*cos(E))-a*e #Asteroid position with x axis corresponding to the perihelion
    r.y = a*sqrt(1-(e**2))*sin(E)
    r.z = 0
    r = rotate(r, angle=-1*wprime, axis=vector(0,0,1)) #Converts position to ecliptic coordinates
    r = rotate(r, angle= iprime, axis=vector(1,0,0))
    r = rotate(r, angle= Oprime, axis=vector(0,0,1))
    return r
    
#Orbital elements of the asteroid
a = 2.1525535780792735
e = 0.455831873080013
M = 0.10125832862508845
Oprime = 3.3741921481515096
iprime = 1.9203619012139899
wprime = 6.27828775408915

#Orbital elements of the Earth
aEarth = 1.00000011  
eEarth = 0.01671022   
IEarth = radians(0.00005) 
OmegaEarth = radians(-11.26064)
wEarth = radians(102.94719) 
MEarth = radians(100.46435)

#More required parameters
mu = 1.000000
time = 2458306.749988
dt = .05 #Time step of visualization

#Parameters needed for the asteroid
period = 2*pi*sqrt(a**3/mu)
r1ecliptic = vector(0, 0, 0) #Creates a position vector for the asteroid
Mtrue = (2*pi/period)*(time) + M
Etrue = solvekep(Mtrue)
r1ecliptic = locateAsteroid(wprime, Oprime, iprime, a, e, Etrue, r1ecliptic) #Calculates position

#Parameters needed for the Earth
periodEarth = 2*pi*sqrt(aEarth**3/mu)
r1Earth = vector(0, 0, 0) #Creates a position vector for the Earth
MtrueEarth = (2*pi/periodEarth)*(time) + M
EtrueEarth = solvekep(MtrueEarth)
r1Earth = locateAsteroid(wEarth, OmegaEarth, IEarth, aEarth, eEarth, EtrueEarth, r1Earth) #Calculates position

# Create objects 
asteroid = sphere(pos=r1ecliptic*150, radius=(10), color=color.white) #Sets up orbit visualization
asteroid.trail = curve(color=color.white)
earth = sphere(pos=r1Earth*150, radius=(15), color=color.blue) #Sets up orbit visualization
earth.trail = curve(color=color.blue)
sun = sphere(pos=vector(0,0,0), radius=(50), color=color.yellow)

# Object labels
asteroidLabel = label(pos = asteroid.pos,
    text='Asteroid', xoffset=20,
    yoffset=20, space=30,
    height=16, border=4,
    font='sans')
earthLabel = label(pos = earth.pos,
    text='Earth', xoffset=20,
    yoffset=20, space=30,
    height=16, border=4,
    font='sans')
sunLabel = label(pos = sun.pos,
    text='Sun', xoffset=20,
    yoffset=20, space=30,
    height=16, border=4,
    font='sans')

while True: #Updates position of asteroid and the Earth
    rate(200) 
    time = time + dt
    
    Mtrue = (2*pi/period)*(time) + M
    Etrue = solvekep(Mtrue)
    r1ecliptic = locateAsteroid(wprime, Oprime, iprime, a, e, Etrue, r1ecliptic)
    asteroid.pos = r1ecliptic*150
    asteroid.trail.append(pos=asteroid.pos)
    asteroidLabel.pos = asteroid.pos
    
    MtrueEarth = (2*pi/periodEarth)*(time) + M
    EtrueEarth = solvekep(MtrueEarth)
    r1Earth = locateAsteroid(wEarth, OmegaEarth, IEarth, aEarth, eEarth, EtrueEarth, r1Earth)
    earth.pos = r1Earth*150
    earth.trail.append(pos=earth.pos)
    earthLabel.pos = earth.pos
