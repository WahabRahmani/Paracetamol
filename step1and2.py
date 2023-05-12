import numpy as np

stream = [100, 0, 0, 0]                         # first two values are useful products, 3rd value is the uknown chemical that is separated in a later step, 4th value is all waste
SUM = sum(stream)
totalcost = 0
totalenvironmentcost = 0
tempcost = 10
volumecost = 1/400


#Modeling reactors

######################################################################
#Reactor 1 (step1)


#Ea and K

step1solvent = [1, 2, 3, 4]
Ea11 = [9000, 9100, 8800, 9200]           #j/mol
Ea12 = [9000, 9300, 8400, 8800]            #j/mol
step1cost = [100, 100, 250, 90]                 #£/flowrate
step1environmentcost = [100, 90, 200, 65]       #/flowrate

K11 = 0 
K12 = 0

Ea11[0] = 9000                 #j/mol
Ea12[0] = 9000                 #j/mol
FA0 = 0.001*0.1*0.5/60          #mol/s           pinene flowrate
T = 293                         #K

X = 0.99
totalflow = 0.001*30.1/60       #l/sec
CA0 = FA0/totalflow             #mol/l
V = 9*0.001                     #l
tau = V/totalflow               #sec            volume/total flowrate


K11 = -(np.log(1-X))/(tau*(np.exp(-Ea11[0]/(8.314*T))+(0.15/0.85)*np.exp(-Ea12[0]/(8.314*T))))
K12 = 0.15*K11/0.85
#print(K11, K12)


def step1(solvent, stream, T, X=0.99):
    #print("STEP 1")

    Ea1 = Ea11[solvent-1]
    Ea2 = Ea12[solvent-1]
    cost = step1cost[solvent-1]
    environmentcost = step1environmentcost[solvent-1]

    V = (-stream[0]*np.log(1-X))/(CA0*(K11*np.exp(-Ea1/(8.314*T))+K12*np.exp(-Ea2/(8.314*T))))
    #print("V=", V)
    Y = (K11*np.exp(-Ea1/(8.314*T)))/(K11*np.exp(-Ea1/(8.314*T))+K12*np.exp(-Ea2/(8.314*T)))
    #print("Y=", Y)

    global totalcost, totalenvironmentcost
    step1totalcost = sum(stream)*cost + tempcost*T + volumecost*V
    step1totalenvironmentcost = sum(stream)*environmentcost
    
    totalcost += step1totalcost
    totalenvironmentcost += step1totalenvironmentcost

    stream[3] = stream[0]*(1-X) + stream[0]*X*(1-Y)
    stream[0] = stream[0]*X*Y
    #print(stream)
    #print("cost: ", step1totalcost)
    #print("Env cost: ", step1totalenvironmentcost)

def step1sep(stream, R):
    stream[3] = 0
    stream[0] = stream[0]*R

#step1(step1solvent[0], stream, 500)
#step1(step1solvent[1], stream, T)
#step1(step1solvent[2], stream, T)
#step1sep(stream, 0.962)

##################################################################

#Reactor 2 (step2)

step2solvent = [1, 2, 3, 4]
Ea21 = [5000, 4500, 4000, 6000]
Ea22 = [5000, 4300, 5000, 6500]
Ea23 = [5000, 5700, 5500, 6500]
step2cost = [100, 200, 120, 60]
step2environmentcost = [100, 150, 90, 60]

K21 = 0
K22 = 0
K23 = 0

totalflow2 = 0.11*0.001/60
FA02 = 6.8*10**-7
T2 = 318
CA02 = FA02/totalflow2
X2 = 0.99

tau2 = 120*60

K21 = -(np.log(1-X2))/(tau2*(np.exp(-Ea21[0]/(8.314*T2))+(0.2/0.6)*np.exp(-Ea22[0]/(8.314*T2))+(0.2/0.6)*np.exp(-Ea23[0]/(8.314*T2))))
K22 = 0.2*K21/0.6
K23 = 0.2*K21/0.6

#print(K21, K22, K23)


def step2(solvent, stream, T, X=0.99):
    #print("STEP 2")

    Ea1 = Ea21[solvent-1]
    Ea2 = Ea22[solvent-1]
    Ea3 = Ea23[solvent-1]
    cost = step2cost[solvent-1]
    environmentcost = step2environmentcost[solvent-1]

    V = (-stream[0]*np.log(1-X2))/(CA02*(K21*np.exp(-Ea1/(8.314*T))+K22*np.exp(-Ea2/(8.314*T))+K23*np.exp(-Ea3/(8.314*T))))
    #print("V=", V)

    Y1 = K21*np.exp(-Ea1/(8.314*T))/(K21*np.exp(-Ea1/(8.314*T))+K22*np.exp(-Ea2/(8.314*T))+K23*np.exp(-Ea3/(8.314*T)))
    Y2 = K22*np.exp(-Ea2/(8.314*T))/(K21*np.exp(-Ea1/(8.314*T))+K22*np.exp(-Ea2/(8.314*T))+K23*np.exp(-Ea3/(8.314*T)))
    #print("Y=", Y1, Y2)

    global totalcost, totalenvironmentcost
    step2totalcost = sum(stream)*cost + tempcost*T + volumecost*V
    step2totalenvironmentcost = sum(stream)*environmentcost

    totalcost += step2totalcost
    totalenvironmentcost += step2totalenvironmentcost

    stream[3] = stream[0]*(1-X)
    stream[2] = stream[0]*X*(1-Y1-Y2)
    stream[1] = stream[0]*X*Y2
    stream[0] = stream[0]*X*Y1

    #print(stream)
    #print("cost: ", step2totalcost)
    #print("Env cost: ", step2totalenvironmentcost)

def step2sep(stream, R):
    stream[3] = 0
    stream[2] = stream[2]*R
    stream[1] = stream[1]*R
    stream[0] = stream[0]*R

#step2(step2solvent[0], stream, T2)
#step2(step2solvent[1], stream, T2)
#step2(step2solvent[2], stream, T2)
#step2sep(stream, 0.95)

#####################################################################################

def s12(step1solvent, T1, step2solvent, T2):
    global totalcost, totalenvironmentcost
    stream = [100, 0, 0, 0]
    step1(step1solvent, stream, T1)
    step1sep(stream, 0.962)

    step2(step2solvent, stream, T2)
    step2sep(stream, 0.95)

    stream = [100 * 100/(stream[1]+stream[0]), 0, 0, 0]
    SUM = sum(stream)
    totalcost = 0
    totalenvironmentcost = 0

    step1(step1solvent, stream, T1)
    step1sep(stream, 0.962)

    step2(step2solvent, stream, T2)
    step2sep(stream, 0.95)

    YIELD = (stream[0]+stream[1])/(SUM*0.60909)                     #0.60909
    TOTAL_COST = totalcost/46701.7                                  #46701.7
    TOTAL_ENVIRONMENTAL_COST = totalenvironmentcost/29708.9           #29708.9

    #weights
    a = 2
    b = -0.6
    c = -0.4

    objective = a*YIELD + b*TOTAL_COST + c*TOTAL_ENVIRONMENTAL_COST

    #print(YIELD, TOTAL_COST, TOTAL_ENVIRONMENTAL_COST, objective)
    return [YIELD, TOTAL_COST, TOTAL_ENVIRONMENTAL_COST, objective]


#s12(1, 293, 1, 318)
#s12(step1solvent[0], 288, step2solvent[1], 298)

#s12(2, 350,   3, 310)
#s12( 2, 305,   3, 330)
