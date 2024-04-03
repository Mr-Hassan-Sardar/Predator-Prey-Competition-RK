# predatorPrey.py
# Program to run simulation of the predator-prey relationship

def PredatorPrey_RK2(DT=0.001, simLength=12):
    numIterations = int(simLength / DT) + 1
    t = 0

    predator_population = 15
    predator_birth_fraction = 0.01
    predator_death_proportionality_constant = 1.06
    prey_population = 100
    prey_birth_fraction = 2
    prey_death_proportionality_constant = 0.02

    tLst = [t]
    predatorLst = [predator_population]
    preyLst = [prey_population]
    for i in range(1, numIterations):
        t = i * DT
        k1_prey = (prey_birth_fraction * prey_population - 
                   prey_death_proportionality_constant * predator_population * prey_population) * DT
        k1_predator = (predator_birth_fraction * prey_population * predator_population - 
                       predator_death_proportionality_constant * predator_population) * DT

        k2_prey = (prey_birth_fraction * (prey_population + 0.5 * k1_prey) - 
                   prey_death_proportionality_constant * (predator_population + 0.5 * k1_predator) * 
                   (prey_population + 0.5 * k1_prey)) * DT
        k2_predator = (predator_birth_fraction * (prey_population + 0.5 * k1_prey) * 
                       (predator_population + 0.5 * k1_predator) - 
                       predator_death_proportionality_constant * (predator_population + 0.5 * k1_predator)) * DT

        prey_population += k2_prey
        predator_population += k2_predator

        tLst.append(t)
        predatorLst.append(predator_population)
        preyLst.append(prey_population)

    return tLst, predatorLst, preyLst


def PredatorPrey_RK4(DT=0.001, simLength=12):
    numIterations = int(simLength / DT) + 1
    t = 0

    predator_population = 15
    predator_birth_fraction = 0.01
    predator_death_proportionality_constant = 1.06
    prey_population = 100
    prey_birth_fraction = 2
    prey_death_proportionality_constant = 0.02

    tLst = [t]
    predatorLst = [predator_population]
    preyLst = [prey_population]
    for i in range(1, numIterations):
        t = i * DT
        k1_prey = (prey_birth_fraction * prey_population - 
                   prey_death_proportionality_constant * predator_population * prey_population) * DT
        k1_predator = (predator_birth_fraction * prey_population * predator_population - 
                       predator_death_proportionality_constant * predator_population) * DT

        k2_prey = (prey_birth_fraction * (prey_population + 0.5 * k1_prey) - 
                   prey_death_proportionality_constant * (predator_population + 0.5 * k1_predator) * 
                   (prey_population + 0.5 * k1_prey)) * DT
        k2_predator = (predator_birth_fraction * (prey_population + 0.5 * k1_prey) * 
                       (predator_population + 0.5 * k1_predator) - 
                       predator_death_proportionality_constant * (predator_population + 0.5 * k1_predator)) * DT

        k3_prey = (prey_birth_fraction * (prey_population + 0.5 * k2_prey) - 
                   prey_death_proportionality_constant * (predator_population + 0.5 * k2_predator) * 
                   (prey_population + 0.5 * k2_prey)) * DT
        k3_predator = (predator_birth_fraction * (prey_population + 0.5 * k2_prey) * 
                       (predator_population + 0.5 * k2_predator) - 
                       predator_death_proportionality_constant * (predator_population + 0.5 * k2_predator)) * DT

        k4_prey = (prey_birth_fraction * (prey_population + k3_prey) - 
                   prey_death_proportionality_constant * (predator_population + k3_predator) * 
                   (prey_population + k3_prey)) * DT
        k4_predator = (predator_birth_fraction * (prey_population + k3_prey) * 
                       (predator_population + k3_predator) - 
                       predator_death_proportionality_constant * (predator_population + k3_predator)) * DT

        prey_population += (1 / 6) * (k1_prey + 2 * k2_prey + 2 * k3_prey + k4_prey)
        predator_population += (1 / 6) * (k1_predator + 2 * k2_predator + 2 * k3_predator + k4_predator)

        tLst.append(t)
        predatorLst.append(predator_population)
        preyLst.append(prey_population)

    return tLst, predatorLst, preyLst
