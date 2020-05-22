#Simulatneous optimisation and simulation script for bioethanol fermentation through batch and fed-batch model. The optimisation it is performed
# in order to find the optimal flow rate to the fed-batch process for achieving the highest glucose and xylose simulatneous uptake rate


#import the packages
from scipy.integrate import odeint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#import the parameters of the variables uptake, inhibitory constants and yields
df = pd.read_excel("Kinetics_overview.xlsx")

#The script it consist is one class (the model itsefl) with 6 functions plus the definition of initial conditions and parameters
class model():
    def __init__(self):
        #In here all the parameters, stoichiometric matrix, initial conditions, process vector ,rates vector, time vector and flow vector
        #are going to be either initialized or defined


        #Parameters

        self.nuMaxGlu = float(df['Value2'][df.index[df['Parameters'] == "nuMaxGlu"]])  # s-1
        self.nuMaxXyl = float(df['Value2'][df.index[df['Parameters'] == "nuMaxXyl"]])  # s-1
        self.Ks_Glu = float(df['Value2'][df.index[df['Parameters'] == "Ks_Glu"]])  # kg Glu m-3
        self.Ks_Xyl = float(df['Value2'][df.index[df['Parameters'] == "Ks_Xyl"]])  # kg Xyl m-3
        self.Ki_Glu = float(df['Value2'][df.index[df['Parameters'] == "Ki_Glu"]])  # kg Glu m-3
        self.Ki_Xyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_Xyl"]])  # kg Xyl m-3
        self.Ki_GluXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_GluXyl"]])  # kg Glu m-3
        self.Y_XGlu = float(df['Value2'][df.index[df['Parameters'] == "Y_XGlu"]])  # kg X/kg Glu
        self.Y_XXyl = float(df['Value2'][df.index[df['Parameters'] == "Y_XXyl"]])  # kg X/kg Xyl
        self.Ki_EtOHmaxGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_EtOHmaxGlu"]])  # kg Glu m-3
        self.Ki_EtOHmaxXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_EtOHmaxXyl"]])  # kg Xyl m-3
        self.Y_EtOHGlu = float(df['Value2'][df.index[df['Parameters'] == "Y_EtOHGlu"]])  # kg EtOH/kg Glu
        self.Y_EtOHXyl = float(df['Value2'][df.index[df['Parameters'] == "Y_EtOHXyl"]])  # kg EtOH/kg Xyl
        self.gammaG = float(df['Value2'][df.index[df['Parameters'] == "gammaG"]])  # no unit
        self.gammaX = float(df['Value2'][df.index[df['Parameters'] == "gammaX"]])  # no unit

        # Acetate parameters
        self.nuHAcMax = float(df['Value2'][df.index[df['Parameters'] == "nuHAcMax"]])  # s-1
        self.Ks_HAc = float(df['Value2'][df.index[df['Parameters'] == "Ks_HAc"]])  # kg HAc m-3
        self.Ki_HAcGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_HAcGlu"]])  # kg HAc m-3
        self.Ki_HAcXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_HAcXyl"]])  # kg HAc m-3
        self.Y_HAcHMF = float(df['Value2'][df.index[df['Parameters'] == "Y_HAcHMF"]])  # kg Ac/kg HMF

        # Furfural parameters
        self.nuFurMax = float(df['Value2'][df.index[df['Parameters'] == "nuFurMax"]])  # s-1
        self.Ks_Fur = float(df['Value2'][df.index[df['Parameters'] == "Ks_Fur"]])  # kg Furfural m-3
        self.Ki_FurGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_FurGlu"]])  # kg Furfural m-3
        self.Ki_FurXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_FurXyl"]])  # kg Furfural m-3
        self.Ki_FurHMF = float(df['Value2'][df.index[df['Parameters'] == "Ki_FurHMF"]])  # kg Furfural m-3
        self.Y_FurFA = float(df['Value2'][df.index[df['Parameters'] == "Y_FurFA"]])  # kg FA/kg Fur

        # Furfuryl alcohol parameters
        self.Ki_FAGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_FAGlu"]])  # kg FA m-3
        self.Ki_FAXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_FAXyl"]])  # kg FA m-3

        # HMF parameters
        self.nuHMFMax = float(df['Value2'][df.index[df['Parameters'] == "nuHMFMax"]])  # s-1
        self.Ks_HMF = float(df['Value2'][df.index[df['Parameters'] == "Ks_HMF"]])  # kg HMF m-3
        self.Ki_HMFGlu = float(df['Value2'][df.index[df['Parameters'] == "Ki_HMFGlu"]])  # kg HMF m-3
        self.Ki_HMFXyl = float(df['Value2'][df.index[df['Parameters'] == "Ki_HMFXyl"]])  # kg HMF m-3


        # Initialization of the process vector used in the rxn function
        self.rho = np.zeros(5)

        # Initialization of the overall rates vector
        self.r = np.zeros(8)

        #Initialization of stoichiometric matric
        self.s = np.zeros((5, 8))

        self.s[0, 0] = self.Y_XGlu
        self.s[1, 0] = self.Y_XXyl
        self.s[2, 0] = 0
        self.s[3, 0] = 0
        self.s[4, 0] = 0

        self.s[0, 1] = -1
        self.s[1, 1] = 0
        self.s[2, 1] = 0
        self.s[3, 1] = 0
        self.s[4, 1] = 0

        self.s[0, 2] = 0
        self.s[1, 2] = -1
        self.s[2, 2] = 0
        self.s[3, 2] = 0
        self.s[4, 2] = 0

        self.s[0, 3] = self.Y_EtOHGlu
        self.s[1, 3] = self.Y_EtOHXyl
        self.s[2, 3] = 0
        self.s[3, 3] = 0
        self.s[4, 3] = 0

        self.s[0, 4] = 0
        self.s[1, 4] = 0
        self.s[2, 4] = -1
        self.s[3, 4] = 0
        self.s[4, 4] = 0

        self.s[0, 5] = 0
        self.s[1, 5] = 0
        self.s[2, 5] = 0
        self.s[3, 5] = -1
        self.s[4, 5] = self.Y_HAcHMF

        self.s[0, 6] = 0
        self.s[1, 6] = 0
        self.s[2, 6] = 0
        self.s[3, 6] = 0
        self.s[4, 6] = -1

        self.s[0, 7] = 0
        self.s[1, 7] = 0
        self.s[2, 7] = self.Y_FurFA
        self.s[3, 7] = 0
        self.s[4, 7] = 0

        #Definition of the initial conditions of the process to be optimized (V0 is going to be the volume considered as initial volume of fedbatch process, then
        # the batch volume is going to be calculated as the last volume of fed-batch process in order to have comparable results
        self.V0 = 100  # L
        self.X_mass = 100  # g
        self.Glu0 = 40  # g/L
        self.Xyl0 = 20  # g/L
        self.EtOH0 = 0  # g/L
        self.Fur0 = 1  # g/L
        self.HAc0 = 1  # g/L
        self.HMF0 = 0.5  # g/L
        self.FA0 = 0  # g/L
        self.Max_Ferm_time = 50  # h
        self.F_max = 150  # maximum operative flow
        self.n = 151  # n of different flows in the flow vector

        # Time vector to solve the ODE. It goes from 0 to the fermentation time as we are
        # looking for one hour intervals to choose the optimum F for every hour
        self.t = np.linspace(0, self.Max_Ferm_time, self.Max_Ferm_time + 1)
        self.t = np.around(self.t, decimals=2)

        #Flow vector as options of control flows to apply to the fed-batch inlet (used for the optimization
        self.F = np.linspace(0, self.F_max, self.n)  # L/h

        # _______________________________________________________________________________________________________________________

#PART 1: Optimisation of the fed-batch process

    #rxn_opti funciton: thi function is used by the odeint solving method in function opti(). In here, the model equations
    # such as uptake equations and mass balances are writen and the return of the function will be the differential equation of the variables respect the time.
    #It is calleed rxn because mainly it is multiplied the stoichiometric matrix by the process vector gicing the r (uptake) and then the mass balance can be performed.
    #Note: This function it is specific for optimisation as the differential equations are used with a loop in where every different flow is going to be tried.
    def rxn_opti(self, C, t):

        X, Glu, Xyl, EtOH, Fur, HAc, HMF, FA, V = C

        # Glucose uptake process
        self.rho[0] = self.nuMaxGlu * X * (Glu / (self.Ks_Glu + Glu + ((Glu ** 2) / self.Ki_Glu)) *
                                           (1 - (EtOH / self.Ki_EtOHmaxGlu) ** self.gammaG) *
                                           (1 / (1 + (Fur / self.Ki_FurGlu))) *
                                           (1 / (1 + (HAc / self.Ki_HAcGlu))) *
                                           (1 / (1 + (HMF / self.Ki_HMFGlu))) *
                                           (1 / (1 + (FA / self.Ki_FAGlu))))
        # Xylose uptake process
        self.rho[1] = self.nuMaxXyl * X * (Xyl / (self.Ks_Xyl + Xyl + ((Xyl ** 2) / self.Ki_Xyl)) *
                                           (1 - (EtOH / self.Ki_EtOHmaxXyl) ** self.gammaX) *
                                           (1 / (1 + (Fur / self.Ki_FurXyl))) *
                                           (1 / (1 + (HAc / self.Ki_HAcXyl))) *
                                           (1 / (1 + (HMF / self.Ki_HMFXyl))) *
                                           (1 / (1 + (FA / self.Ki_FAXyl))) *
                                           (1 / (1 + (Glu / self.Ki_GluXyl))))
        # Fur uptake process
        self.rho[2] = self.nuFurMax * X * (Fur / (self.Ks_Fur + Fur))
        # HAc uptake process
        self.rho[3] = self.nuHAcMax * X * (HAc / (self.Ks_HAc + HAc))
        # HMF uptake process
        self.rho[4] = self.nuHMFMax * X * (HMF / (self.Ks_HMF + HMF)) * (1 / (1 + (Fur / self.Ki_FurGlu)))

        # Overall conversion rate (stoichiometric matrix * process vector)
        self.r[0] = self.s[0, 0] * self.rho[0] + self.s[1, 0] * self.rho[1] + self.s[2, 0] * self.rho[2] + \
                    self.s[
                        3, 0] * \
                    self.rho[3] + self.s[4, 0] * self.rho[4]
        self.r[1] = self.s[0, 1] * self.rho[0] + self.s[1, 1] * self.rho[1] + self.s[2, 1] * self.rho[2] + \
                    self.s[
                        3, 1] * \
                    self.rho[3] + self.s[4, 1] * self.rho[4]
        self.r[2] = self.s[0, 2] * self.rho[0] + self.s[1, 2] * self.rho[1] + self.s[2, 2] * self.rho[2] + \
                    self.s[
                        3, 2] * \
                    self.rho[3] + self.s[4, 2] * self.rho[4]
        self.r[3] = self.s[0, 3] * self.rho[0] + self.s[1, 3] * self.rho[1] + self.s[2, 3] * self.rho[2] + \
                    self.s[
                        3, 3] * \
                    self.rho[3] + self.s[4, 3] * self.rho[4]
        self.r[4] = self.s[0, 4] * self.rho[0] + self.s[1, 4] * self.rho[1] + self.s[2, 4] * self.rho[2] + \
                    self.s[
                        3, 4] * \
                    self.rho[3] + self.s[4, 4] * self.rho[4]
        self.r[5] = self.s[0, 5] * self.rho[0] + self.s[1, 5] * self.rho[1] + self.s[2, 5] * self.rho[2] + \
                    self.s[
                        3, 5] * \
                    self.rho[3] + self.s[4, 5] * self.rho[4]
        self.r[6] = self.s[0, 6] * self.rho[0] + self.s[1, 6] * self.rho[1] + self.s[2, 6] * self.rho[2] + \
                    self.s[
                        3, 6] * \
                    self.rho[3] + self.s[4, 6] * self.rho[4]
        self.r[7] = self.s[0, 7] * self.rho[0] + self.s[1, 7] * self.rho[1] + self.s[2, 7] * self.rho[2] + \
                    self.s[
                        3, 7] * \
                    self.rho[3] + self.s[4, 7] * self.rho[4]

        # Mass balances

        dVdt = self.F[self.i]

        dXdt = self.r[0] - X * (self.F[self.i] / V)

        dGludt = self.r[1] + (self.F[self.i] / V) * (self.Glu0 - Glu)

        dXyldt = self.r[2] + (self.F[self.i] / V) * (self.Xyl0 - Xyl)

        dEtOHdt = self.r[3] - EtOH * (self.F[self.i] / V)

        dFurdt = self.r[4] + (self.F[self.i] / V) * (self.Fur0 - Fur)

        dHAcdt = self.r[5] + (self.F[self.i] / V) * (self.HAc0 - HAc)

        dHMFdt = self.r[6] + (self.F[self.i] / V) * (self.Fur0 - Fur)

        dFAdt = self.r[7] - FA * (self.F[self.i] / V)

        return [dXdt, dGludt, dXyldt, dEtOHdt, dFurdt, dHAcdt, dHMFdt, dFAdt,
                dVdt]  # Here, we are solving the differential equation for each point.


    #opti funciton:
    #this function is the one required to solve the optimisation. It uses the odeint method from scipy to get
    #solve the differential equations of the value of the state variables at every defined point.
    #The optimisation works with two loops: the first one loops all the function for every hour, so if the max time fermentation
    #has been set at 50 ours it will perform 50 loops. For each of this first loop it will try every different flow in the flow
    #defined in the beginning. After every hour and all the flows already tried, the system will save as the initial condition
    #for next hour simulation the values of the state variables that makes the maximum value of the
    #objective function (glucose + xylose uptake rate). The function will save the Final volume of the fed-batch
    # (so can be used for batch comparison) and the optimal flow profile to be applied in the fed-batch solver
    # (which is printed in a plot)

    def opti(self):
        # The C0 matrix is created in order to update each of the rows with the values of the simulation which gives the optimum Q
        Cinit = np.zeros((len(self.t), 9))

        # Vector with the initial conditions of the state variables


        self.X0 = self.X_mass / self.V0  # g/L

        Cinit[0] = [self.X0, self.Glu0, self.Xyl0, self.EtOH0, self.Fur0, self.HAc0, self.HMF0, self.FA0, self.V0]

        # Vector that will save the values of Xylose rates  at each Q and each timepoint evaluation

        rates = np.zeros((len(self.t), len(self.F)))

        # This Vector will be the vector that contains the optimal Q for every timepoint
        self.F_opt_rates = np.zeros(len(self.t))

        #initialize the variable Vfinal which will set the total volume of fermentation
        self.Vfinal = 0

        for j in range(len(self.t)):
            # Then the next "for" is created to try all the different proposed Q from the Q vectors

            for i in range(len(self.F)):
                self.i = i
                # Here is performed the solving method: first a vector t is created in order to perform de simulation at every timpoint (j)

                timepoints = np.linspace(j, j + 1, 2)
                timepoints = np.around(timepoints, decimals=2)

                # C vector containing the result of the solve ODEint
                C = odeint(self.rxn_opti, Cinit[j], timepoints)
                # Now it is going to tart the collect of optimal F
                # First there is a condition where if the timepoint arrives to the last time point the script breaks just for a matter of how the next steps are defined
                # (In case to include this timepoint Ferm_time should be one more than actual value)

                if j == (self.Max_Ferm_time):
                    break
                # For the first timepoint the initial Xyl rate uses the first value of Xyl -->Xyl0
                # The xylose rates are going to be measured with he last values of the simulation "C[-1,...]" following the mass balance
                # which also takes into account the dilution of the added F

                rates[j, i] = self.Y_EtOHXyl * (
                        self.nuMaxXyl * C[-1, 0] * (
                            C[-1, 2] / (self.Ks_Xyl + C[-1, 2] + ((C[-1, 2] ** 2) / self.Ki_Xyl)) *
                            (1 - (C[-1, 3] / self.Ki_EtOHmaxXyl) ** self.gammaX) *
                            (1 / (1 + (C[-1, 4] / self.Ki_FurXyl))) *
                            (1 / (1 + (C[-1, 5] / self.Ki_HAcXyl))) *
                            (1 / (1 + (C[-1, 6] / self.Ki_HMFXyl))) *
                            (1 / (1 + (C[-1, 7] / self.Ki_FAXyl))) *
                            (1 / (1 + (C[-1, 1] / self.Ki_GluXyl))))) + \
                              self.Y_EtOHGlu * (self.nuMaxGlu * C[-1, 0] * (
                        C[-1, 1] / (self.Ks_Glu + C[-1, 1] + ((C[-1, 2] ** 1) / self.Ki_Glu)) *
                        (1 - (C[-1, 3] / self.Ki_EtOHmaxGlu) ** self.gammaG) *
                        (1 / (1 + (C[-1, 4] / self.Ki_FurGlu))) *
                        (1 / (1 + (C[-1, 5] / self.Ki_HAcGlu))) *
                        (1 / (1 + (C[-1, 6] / self.Ki_HMFGlu))) *
                        (1 / (1 + (C[-1, 7] / self.Ki_FAGlu)))))

                # #A dataframe is created for an easier visualization
                # Here comes the main point of the optimization. the initial concentrations of the next timepoint will be the values achieved with the simulation of the optimal F
                if rates[j, i] == rates[j].max():
                    Cinit[j + 1] = C[-1]

                # Here a restriction of volume is included

                # #The array for the optimal F at ach timepoint is created below
                # for m in range(len(timepoints)):

                self.F_opt_rates[j] = self.F[np.argmax(rates, axis=1)[j]]

                if Cinit[j, 8] >= 8.5 * Cinit[0, 8]:
                    self.F_opt_rates[j] = 0
                    self.Vfinal = float(Cinit[np.where(Cinit[:, 8] >= 8.5 * Cinit[0, 8]), 8][:, 0])

        self.F_opt_rates_data = pd.DataFrame(self.F_opt_rates, index=self.t)

        print("Final volume of reaction = " + str(self.Vfinal) + " L")

        # Plotting
        plt.plot(self.t, self.F_opt_rates)
        plt.title("Optimal Flow vs time")
        plt.xlabel("Time (h)")
        plt.ylabel("Flow (L/h)")
        plt.savefig('Optimalflow.png')
        plt.show()

        return
    # ______________________________________________________________________________________________________________________

#PART 2: Simulation

    #Flow function:
    #Function used as an argument for the odeint solver function of the fed-batch process as we will need one specific flow for
    #time solving the equation and this will come from the F_opt_rates_data calculated in the optimisation
    def Flow(self,t):



        return self.F_opt_rates_data[0][round(t)]  # L/h

    #rxn function:
    #Same function as rxn_opti, but in this case, as it is used for the simulation the F is not going to be looped and only
    #the appropiate Flow will be used for the fed-batch solver, and a F=0 will be used for batch solver
    def rxn(self, C, t, F):
        # State variables inside of the Vector C which contains the concentration of the variables at
        # every single time

        X, Glu, Xyl, EtOH, Fur, HAc, HMF, FA, V = C




        if self.Mode == "Batch":
            self.F = 0
        if self.Mode == "FedBatch":
            self.F = self.Flow(t)

        # Glucose uptake process
        self.rho[0] = self.nuMaxGlu * X * (Glu / (self.Ks_Glu + Glu + ((Glu ** 2) / self.Ki_Glu)) *
                                           (1 - (EtOH / self.Ki_EtOHmaxGlu) ** self.gammaG) *
                                           (1 / (1 + (Fur / self.Ki_FurGlu))) *
                                           (1 / (1 + (HAc / self.Ki_HAcGlu))) *
                                           (1 / (1 + (HMF / self.Ki_HMFGlu))) *
                                           (1 / (1 + (FA / self.Ki_FAGlu))))
        # Xylose uptake process
        self.rho[1] = self.nuMaxXyl * X * (Xyl / (self.Ks_Xyl + Xyl + ((Xyl ** 2) / self.Ki_Xyl)) *
                                           (1 - (EtOH / self.Ki_EtOHmaxXyl) ** self.gammaX) *
                                           (1 / (1 + (Fur / self.Ki_FurXyl))) *
                                           (1 / (1 + (HAc / self.Ki_HAcXyl))) *
                                           (1 / (1 + (HMF / self.Ki_HMFXyl))) *
                                           (1 / (1 + (FA / self.Ki_FAXyl))) *
                                           (1 / (1 + (Glu / self.Ki_GluXyl))))
        # Fur uptake process
        self.rho[2] = self.nuFurMax * X * (Fur / (self.Ks_Fur + Fur))
        # HAc uptake process
        self.rho[3] = self.nuHAcMax * X * (HAc / (self.Ks_HAc + HAc))
        # HMF uptake process
        self.rho[4] = self.nuHMFMax * X * (HMF / (self.Ks_HMF + HMF)) * (1 / (1 + (Fur / self.Ki_FurGlu)))

        # Overall conversion rate (stoichiometric matrix * process vector)
        self.r[0] = self.s[0, 0] * self.rho[0] + self.s[1, 0] * self.rho[1] + self.s[2, 0] * self.rho[2] + self.s[
            3, 0] * \
                    self.rho[3] + self.s[4, 0] * self.rho[4]
        self.r[1] = self.s[0, 1] * self.rho[0] + self.s[1, 1] * self.rho[1] + self.s[2, 1] * self.rho[2] + self.s[
            3, 1] * \
                    self.rho[3] + self.s[4, 1] * self.rho[4]
        self.r[2] = self.s[0, 2] * self.rho[0] + self.s[1, 2] * self.rho[1] + self.s[2, 2] * self.rho[2] + self.s[
            3, 2] * \
                    self.rho[3] + self.s[4, 2] * self.rho[4]
        self.r[3] = self.s[0, 3] * self.rho[0] + self.s[1, 3] * self.rho[1] + self.s[2, 3] * self.rho[2] + self.s[
            3, 3] * \
                    self.rho[3] + self.s[4, 3] * self.rho[4]
        self.r[4] = self.s[0, 4] * self.rho[0] + self.s[1, 4] * self.rho[1] + self.s[2, 4] * self.rho[2] + self.s[
            3, 4] * \
                    self.rho[3] + self.s[4, 4] * self.rho[4]
        self.r[5] = self.s[0, 5] * self.rho[0] + self.s[1, 5] * self.rho[1] + self.s[2, 5] * self.rho[2] + self.s[
            3, 5] * \
                    self.rho[3] + self.s[4, 5] * self.rho[4]
        self.r[6] = self.s[0, 6] * self.rho[0] + self.s[1, 6] * self.rho[1] + self.s[2, 6] * self.rho[2] + self.s[
            3, 6] * \
                    self.rho[3] + self.s[4, 6] * self.rho[4]
        self.r[7] = self.s[0, 7] * self.rho[0] + self.s[1, 7] * self.rho[1] + self.s[2, 7] * self.rho[2] + self.s[
            3, 7] * \
                    self.rho[3] + self.s[4, 7] * self.rho[4]

        # Mass balances

        dVdt = self.F

        dXdt = self.r[0] - X * (self.F / V)

        dGludt = self.r[1] + (self.F / V) * (self.Glu0 - Glu)

        dXyldt = self.r[2] + (self.F / V) * (self.Xyl0 - Xyl)

        dEtOHdt = self.r[3] - EtOH * (self.F / V)

        dFurdt = self.r[4] + (self.F / V) * (self.Fur0 - Fur)

        dHAcdt = self.r[5] + (self.F / V) * (self.HAc0 - HAc)

        dHMFdt = self.r[6] + (self.F / V) * (self.Fur0 - Fur)

        dFAdt = self.r[7] - FA * (self.F / V)


        return [dXdt, dGludt, dXyldt, dEtOHdt, dFurdt, dHAcdt, dHMFdt, dFAdt,
                dVdt]  # Here, we are solving the differential equation for each point.

    #solve_batch function:
    #simulates the batch process using as volume the final volume of the optimisation to be comparable to the fed-batch.
    #Flow is set to 0 with the variable Mode=Batch. It returns the plots of the simulation
    def solve_batch(self):


        V0 = float(self.Vfinal)  # L
        self.X0 = self.X_mass / V0  # g/L
        self.C0 = [self.X0, self.Glu0, self.Xyl0, self.EtOH0, self.Fur0, self.HAc0, self.HMF0, self.FA0, V0]


        self.Mode = "Batch"

        # C vector containing the result of the solve
        C = odeint(self.rxn, self.C0, self.t, args=(self.Flow,))
        t_fin = np.where(C[:, 2] < 0.5)[0]

        print("Batch Volume: " + str(V0) + " L")
        print("Final time of fermentation batch mode: " + str(float(self.t[t_fin[0]])) + " hours")


        # plotting code
        plt.plot(self.t, C[:, 0], 'b-', label="X")
        plt.plot(self.t, C[:, 1], 'g-', label="Glu")
        plt.plot(self.t, C[:, 2], 'r-', label="Xyl")
        plt.plot(self.t, C[:, 3], 'c-', label="EtOH")
        plt.plot(self.t, C[:, 4], 'm-', label="Fur")
        plt.plot(self.t, C[:, 5], 'y-', label="HAc")
        plt.plot(self.t, C[:, 6], 'k-', label="HMF")
        plt.plot(self.t, C[:, 7], '0.75', label="FA")
        plt.xlabel("Time (h)")
        plt.ylabel("Conc. (g/L)")
        plt.legend(loc="upper right")
        plt.title('Batch simulation')
        plt.savefig('Sim_batch.png')
        plt.show()

    #solve_fed-batch function:
    #simulates the fedbatch process.
    #It returns the plots of the simulation and a plot showing the optimal flow rate and the volume vs time
    def solve_fedbatch(self):


        self.X0 = self.X_mass / self.V0  # g/L
        self.C0 = [self.X0, self.Glu0, self.Xyl0, self.EtOH0, self.Fur0, self.HAc0, self.HMF0, self.FA0, self.V0]

        self.Mode="FedBatch"
        # Vector containing all the possible values of F to detect the optimal value from 0 to Fmax and n will define the amount of values to evaluate

        # C vector containing the result of the solve
        C = odeint(self.rxn, self.C0, self.t, args=(self.Flow,))
        t_fin = np.where(C[:, 2] < 0.5)[0]

        print("Initial volume Fed-batch Volume: " + str(self.V0) + " L")
        print("Final volume Fed-batch Volume: " + str(float(C[-1,8])) + " L")
        print("Final time of fermentation fed-batch mode : " + str(float(self.t[t_fin[0]])) + " hours")


        # plotting code
        plt.plot(self.t, C[:, 0], 'b-', label="X")
        plt.plot(self.t, C[:, 1], 'g-', label="Glu")
        plt.plot(self.t, C[:, 2], 'r-', label="Xyl")
        plt.plot(self.t, C[:, 3], 'c-', label="EtOH")
        plt.plot(self.t, C[:, 4], 'm-', label="Fur")
        plt.plot(self.t, C[:, 5], 'y-', label="HAc")
        plt.plot(self.t, C[:, 6], 'k-', label="HMF")
        plt.plot(self.t, C[:, 7], '0.75', label="FA")
        plt.xlabel("Time (h)")
        plt.ylabel("Conc. (g/L)")
        plt.legend(loc="upper right")
        plt.title('FedBatch simulation')
        plt.savefig('Sim_fedbatch.png')
        plt.show()

        # plotting to see the volume vs the flow
        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('time (hours)')
        ax1.set_ylabel('F opt for max rates in (L/h)', color=color)
        ax1.plot(self.t, self.F_opt_rates, color=color)
        ax1.tick_params(axis='y', labelcolor=color)

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        color = 'tab:blue'
        ax2.set_ylabel('Volume', color=color)  # we already handled the x-label with ax1
        ax2.plot(self.t, C[:, 8], color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.title("Flow and Volume vs time of fermentation")
        plt.savefig('Flow_volume_sim_fedbatch.png')
        plt.show()


#Here the model is solved:
#for that, first the optimisation has to be initialized in order to simulate the batch and fed-batch process


M = model()
M.opti()
M.solve_batch()
M.solve_fedbatch()



