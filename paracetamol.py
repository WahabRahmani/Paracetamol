#paracetamol

from matplotlib import pyplot as plt
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

    #return V

def step1sep(stream, R):
    stream[3] = 0
    stream[0] = stream[0]*R

#step1(step1solvent[0], [100, 0, 0, 0], 330)
#step1(step1solvent[0], [100, 0, 0, 0], 320)
#step1(step1solvent[0], [100, 0, 0, 0], 300)
#step1(step1solvent[0], stream, 290)
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

#step2(step2solvent[0], stream, 318)
#step2(step2solvent[0], [100, 0, 0, 0], 380)
#step2(step2solvent[0], [100, 0, 0, 0], 300)
#step2(step2solvent[0], [100, 0, 0, 0], 350)
#step2sep(stream, 0.95)

#####################################################################################

#Reactor 3 (step3)


#Ea and K

step3solvent = [1, 2, 3, 4]
Ea31 = [9000, 9100, 7500, 9500]           #j/mol
Ea32 = [9000, 8600, 7900, 10000]            #j/mol
step3cost = [100, 90, 250, 100]                 #£/flowrate
step3environmentcost = [100, 90, 200, 60]       #/flowrate

K31 = 0 
K32 = 0

Ea31[0] = 9000                 #j/mol
Ea32[0] = 9000                 #j/mol
FA03 = 0.001*0.1*0.5/60          #mol/s           pinene flowrate
T3 = 293                         #K

X3 = 0.99
totalflow = 0.001*30.1/60       #l/sec
CA03 = FA03/totalflow             #mol/l
V = 9*0.001                     #l
tau3 = V/totalflow               #sec            volume/total flowrate


K31 = -(np.log(1-X3))/(tau3*(np.exp(-Ea31[0]/(8.314*T3))+(0.15/0.85)*np.exp(-Ea32[0]/(8.314*T3))))
K32 = 0.15*K31/0.85
#print(K31, K32)


def step3(solvent, stream, T, X=0.99):
    #print("STEP 3")

    Ea1 = Ea31[solvent-1]
    Ea2 = Ea32[solvent-1]
    cost = step3cost[solvent-1]
    environmentcost = step3environmentcost[solvent-1]

    V = (-(stream[0]+stream[1])*np.log(1-X))/(CA03*(K31*np.exp(-Ea1/(8.314*T))+K32*np.exp(-Ea2/(8.314*T))))
    #print("V=", V)
    Y = (K31*np.exp(-Ea1/(8.314*T)))/(K31*np.exp(-Ea1/(8.314*T))+K32*np.exp(-Ea2/(8.314*T)))
    #print("Y=", Y)

    global totalcost, totalenvironmentcost
    step3totalcost = sum(stream)*cost + tempcost*T + volumecost*V
    step3total_env_cost = sum(stream)*environmentcost
    
    totalcost += step3totalcost
    totalenvironmentcost += step3total_env_cost

    stream[3] = (stream[0]+stream[1])*(1-X) + (stream[0]+stream[1])*X*(1-Y)
    stream[0] = stream[0]*X*Y
    stream[1] = stream[1]*X*Y

    #print(stream)
    #print("cost: ", step3totalcost)
    #print("Env cost: ", step3total_env_cost)

def step3sep(stream, R):
    stream[3] = 0
    stream[2] = stream[2]*R
    stream[0] = stream[0]*R
    stream[1] = stream[1]*R

#step3(step3solvent[0], [100, 0, 0, 0], 270)
#step3(step3solvent[0], stream, 290)
#step3(step3solvent[0], [100, 0, 0, 0], 320)
#step3(step3solvent[0], [100, 0, 0, 0], 350)
#step3sep(stream, 0.9)

###################################################################################

#Reactor 4 (step4)

step4solvent = [1, 2, 3, 4]
step4catalyst = [1, 2, 3, 4]

Ea41 = [[5000, 4000, 6500, 4200],
        [4800, 3800, 6500, 3700],
        [5300, 4400, 6500, 4400],
        [4900, 3900, 6300, 4500]]

Ea42 = [[5000, 4000, 7000, 5000],
        [4800, 3900, 6700, 4500],
        [5300, 4300, 7000, 4800],
        [4900, 4000, 6500, 4900]]
 
step4solvent_cost = [100, 200, 60, 120]
step4catalyst_cost = [0.0100, 0.0200, 0.0060, 0.0120]
step4solvent_env_cost = [100, 150, 60, 90]
step4catalyst_env_cost = [0.0100, 0.0200, 0.0060, 0.0120]


totalflow4 = (0.29+15.7)*0.001/60           #l/s
tau4 = 12                             
FA04 = 0.001*0.1*0.5/60
T4 = 423
X4 = 0.99
CA04 = FA04/totalflow4

K41 = -(np.log(1-X4))/(tau4*(np.exp(-Ea41[0][0]/(8.314*T4))+(0.05/0.95)*np.exp(-Ea42[0][0]/(8.314*T4))))
K42 = (0.05/0.95)*K41

#print(K41, K42)

def step4(solvent, catalyst, stream, T, X=0.99):
    #print("STEP 4")

    Ea1 = Ea41[catalyst-1][solvent-1]
    Ea2 = Ea42[catalyst-1][solvent-1]

    V = (-(stream[0]+stream[1])*np.log(1-X))/(CA04*(K41*np.exp(-Ea1/(8.314*T))+K42*np.exp(-Ea2/(8.314*T))))
    #print("V=", V)
    Y = (K41*np.exp(-Ea1/(8.314*T)))/(K41*np.exp(-Ea1/(8.314*T))+K42*np.exp(-Ea2/(8.314*T)))
    #print("Y=", Y)

    global totalcost, totalenvironmentcost
    step4totalcost = sum(stream)*step4solvent_cost[solvent-1] + V*step4catalyst_cost[catalyst-1] + tempcost*T + volumecost*V
    step4total_env_cost = sum(stream)*step4solvent_env_cost[solvent-1] + V*step4catalyst_env_cost[catalyst-1]
    
    totalcost += step4totalcost
    totalenvironmentcost += step4total_env_cost

    stream[3] = (stream[0]+stream[1])*(1-X) + (stream[0]+stream[1])*X*(1-Y)
    stream[0] = stream[0]*X*Y
    stream[1] = stream[1]*X*Y

    #print(stream)
    #print("cost: ", step4totalcost)
    #print("Env cost: ", step4total_env_cost)
    #print(V*step4catalyst_cost[catalyst-1])
    #print(sum(stream)*step4solvent_cost[solvent-1])
    #print(tempcost*T)
    #print(volumecost*V)


def step4sep(stream, R1, R2):
    stream[3] = 0
    stream[2] = 0
    stream[0] = stream[0]*R1
    stream[1] = stream[1]*R2

#[34.5967953317415, 11.532265110580498, 13.704414867000002, 0]

#step4(step4solvent[0], step4catalyst[0], stream, T4)
#step4(step4solvent[0], step4catalyst[0], [34.5967953317415, 11.532265110580498, 13.704414867000002, 0], 350)
#step4(step4solvent[0], step4catalyst[0], [34.5967953317415, 11.532265110580498, 13.704414867000002, 0], 400)
#step4(step4solvent[0], step4catalyst[0], [34.5967953317415, 11.532265110580498, 13.704414867000002, 0], 450)

#step4sep(stream, 0.968, 0.979)

###############################################################################

#Reactor 5 (step5)

step5solvent = [1, 2, 3, 4]
step5catalyst = [1, 2, 3, 4]

Ea51 = [[5000, 4000, 6000, 4000],
        [4500, 3600, 5600, 4200],
        [5500, 4000, 5700, 4500],
        [5000, 3900, 5900, 3900]]

Ea52 = [[5000, 4000, 6000, 4500],
        [4700, 3600, 5900, 4700],
        [5300, 4200, 5800, 5300],
        [4700, 3800, 6000, 4300]]

step5solvent_cost = [100, 200, 60, 120]
step5catalyst_cost = [0.100, 0.200, 0.060, 0.120]
step5solvent_env_cost = [100, 150, 60, 90]
step5catalyst_env_cost = [0.0100, 0.0200, 0.0060, 0.0120]

totalflow5 = 20.54*0.001/60           #l/s
tau5 = 60*60                             
FA05 = 0.0658/60
T5 = 383
X5 = 0.98
CA05 = FA05/totalflow5

K51 = -(np.log(1-X5))/(tau5*(np.exp(-Ea51[0][0]/(8.314*T5))+(0.05/0.95)*np.exp(-Ea52[0][0]/(8.314*T5))))
K52 = (0.05/0.95)*K51

#print(K51, K52)

def step5(solvent, catalyst, stream, T, X=0.98):
    #print("STEP 5")

    Ea1 = Ea51[catalyst-1][solvent-1]
    Ea2 = Ea52[catalyst-1][solvent-1]

    V = (-(stream[0])*np.log(1-X))/(CA05*(K51*np.exp(-Ea1/(8.314*T))+K52*np.exp(-Ea2/(8.314*T))))
    #print("V=", V)
    Y = (K51*np.exp(-Ea1/(8.314*T)))/(K51*np.exp(-Ea1/(8.314*T))+K52*np.exp(-Ea2/(8.314*T)))
    #print("Y=", Y)

    global totalcost, totalenvironmentcost
    step5totalcost = sum(stream)*step5solvent_cost[solvent-1] + V*step5catalyst_cost[catalyst-1] + tempcost*T + volumecost*V*10
    step5total_env_cost = sum(stream)*step5solvent_env_cost[solvent-1] + V*step5catalyst_env_cost[catalyst-1]
    
    totalcost += step5totalcost
    totalenvironmentcost += step5total_env_cost

    stream[3] = (stream[0])*(1-X) + (stream[0])*X*(1-Y)
    stream[0] = stream[0]*X*Y

    #print(stream)
    #print("cost: ", step5totalcost)
    #print("Env cost: ", step5total_env_cost)
    #print(V*step5catalyst_cost[catalyst-1])
    #print(sum(stream)*step5solvent_cost[solvent-1])
    #print(tempcost*T)
    #print(volumecost*V*10)

def step5sep(stream, R):
    stream[3] = 0
    stream[2] = 0
    stream[0] = stream[0]*R

#[32.538286009502876, 10.846095336500959, 13.704414867000002, 2.744679096318161]

#step5(step5solvent[0], step5catalyst[0], stream, T5)
#step5(step5solvent[0], step5catalyst[0], [32.538286009502876, 10.846095336500959, 13.704414867000002, 2.744679096318161], 380)
#step5(step5solvent[0], step5catalyst[0], [32.538286009502876, 10.846095336500959, 13.704414867000002, 2.744679096318161], 420)
#step5(step5solvent[0], step5catalyst[0], [32.538286009502876, 10.846095336500959, 13.704414867000002, 2.744679096318161], 460)

#step5sep(stream, 0.931)


###############################################################################

#Reactor 6 (step6)

step6solvent = [1, 2, 3, 4]
step6catalyst = [1, 2, 3, 4]

Ea61 = [[5000, 4000, 6000, 4800],
        [4700, 3900, 5900, 4600],
        [5200, 4100, 6200, 5000],
        [5100, 3900, 5800, 5000]]

Ea62 = [[5000, 4000, 6500, 5000],
        [4800, 3700, 6000, 5200],
        [5100, 4200, 6100, 5400],
        [5000, 4100, 6100, 4900]]

step6solvent_cost = [80, 160, 50, 100]
step6catalyst_cost = [0.0100, 0.0200, 0.0060, 0.0120]
step6solvent_env_cost = [90, 150, 60, 100]
step6catalyst_env_cost = [0.0100, 0.0200, 0.0060, 0.0120]


totalflow6 = 20.54*0.001/60           #l/s
tau6 = 15*60*60                             
FA06 = 0.0658/60
T6 = 503
X6 = 0.93
CA06 = FA06/totalflow6

K61 = -(np.log(1-X6))/(tau6*(np.exp(-Ea61[0][0]/(8.314*T6))+(0.05/0.95)*np.exp(-Ea62[0][0]/(8.314*T6))))
K62 = (0.05/0.95)*K61


def step6(solvent, catalyst, stream, T, X=0.93):
    #print("STEP 6")

    Ea1 = Ea61[catalyst-1][solvent-1]
    Ea2 = Ea62[catalyst-1][solvent-1]

    V = (-(stream[1])*np.log(1-X))/(CA06*(K61*np.exp(-Ea1/(8.314*T))+K62*np.exp(-Ea2/(8.314*T))))
    #print("V=", V)
    Y = (K61*np.exp(-Ea1/(8.314*T)))/(K61*np.exp(-Ea1/(8.314*T))+K62*np.exp(-Ea2/(8.314*T)))
    #print("Y=", Y)

    global totalcost, totalenvironmentcost
    step6totalcost = sum(stream)*step6solvent_cost[solvent-1] + V*step6catalyst_cost[catalyst-1] + tempcost*T/2 + volumecost*V*4
    step6total_env_cost = sum(stream)*step6solvent_env_cost[solvent-1] + V*step6catalyst_env_cost[catalyst-1]
    
    totalcost += step6totalcost
    totalenvironmentcost += step6total_env_cost

    stream[3] = (stream[1])*(1-X) + (stream[1])*X*(1-Y)
    stream[1] = stream[1]*X*Y

    #print(stream)
    #print("cost: ", step6totalcost)
    #print("Env cost: ", step6total_env_cost)
    #print(V*step6catalyst_cost[catalyst-1])
    #print(sum(stream)*step6solvent_cost[solvent-1])
    #print(tempcost*T/2)
    #print(volumecost*V*4)


def step6sep(stream, R):
    stream[3] = 0
    stream[2] = 0
    stream[1] = stream[1]*R


#step6(step6solvent[0], step6catalyst[0], stream, T6)
#step6sep(stream, 0.88)



###############################################################################

# print(stream)
# YIELD = (stream[0]+stream[1])/SUM
# TOTAL_COST = totalcost
# TOTAL_ENVIRONMENTAL_COST = totalenvironmentcost

# print("Total yield:", YIELD)
# print("Total cost:", TOTAL_COST)
# print("Total environmental cost:", TOTAL_ENVIRONMENTAL_COST)

################################################################################


def paracetamol(solvent1, T1, solvent2, T2, solvent3, T3, solvent4, catalyst4, T4, solvent5, catalyst5, T5, solvent6, catalyst6, T6):       #Returns array of yield, cost, env cost, and obj
    global totalcost, totalenvironmentcost
    stream = [100, 0, 0, 0]
    step1(solvent1, stream, T1)
    step1sep(stream, 0.962)

    step2(solvent2, stream, T2)
    step2sep(stream, 0.95)

    step3(solvent3, stream, T3)
    step3sep(stream, 0.9)
    
    step4(solvent4, catalyst4, stream, T4)
    step4sep(stream, 0.968, 0.979)

    step5(solvent5, catalyst5, stream, T5)
    step5sep(stream, 0.931)

    step6(solvent6, catalyst6, stream, T6)
    step6sep(stream, 0.88)

    stream = [100 * (100/sum(stream)), 0, 0, 0]             #To achieve output of 100
    SUM = sum(stream)
    totalcost = 0
    totalenvironmentcost = 0

    step1(solvent1, stream, T1)
    step1sep(stream, 0.962)

    step2(solvent2, stream, T2)
    step2sep(stream, 0.95)

    step3(solvent3, stream, T3)
    step3sep(stream, 0.9)
    
    step4(solvent4, catalyst4, stream, T4)
    step4sep(stream, 0.968, 0.979)

    step5(solvent5, catalyst5, stream, T5)
    step5sep(stream, 0.931)

    step6(solvent6, catalyst6, stream, T6)
    step6sep(stream, 0.88)

    YIELD = (stream[0]+stream[1])/(SUM*0.3555596)
    TOTAL_COST = totalcost/(181127.4)                   #10*181127.4
    TOTAL_ENVIRONMENTAL_COST = totalenvironmentcost/(122650.7 )      #10*122650.7          #RELATIVE to original choice of catalyst and solvent

    #weights
    a = 2
    b = -0.6
    c = -0.4

    objective = a*YIELD + b*TOTAL_COST + c*TOTAL_ENVIRONMENTAL_COST
    #print(YIELD, TOTAL_COST, TOTAL_ENVIRONMENTAL_COST, objective)
    return [YIELD, TOTAL_COST, TOTAL_ENVIRONMENTAL_COST, objective]


#paracetamol(step1solvent[0], T, step2solvent[0], T2, step3solvent[0], T3, step4solvent[0], step4catalyst[0], T4, step5solvent[0], step5catalyst[0], T5, step6solvent[0], step6catalyst[0], T6)
#paracetamol(step1solvent[0], T, step2solvent[0], T2, step3solvent[1], T3, step4solvent[0], step4catalyst[0], T4, step5solvent[0], step5catalyst[0], T5, step6solvent[0], step6catalyst[0], T6)
#paracetamol(2, 330,   3, 310,   4, 288,   3,   3, 430,   3,   3, 400,   3, 3, 540)
#paracetamol(2, 298,   3, 313,   4, 298,   3,   1, 428,   4,   4, 378,   4, 4, 540)
