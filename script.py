#!/usr/bin/env python3
"""
Created on 19/10/2023 [dd/mm/yyyy]
Code by L. J. Duarte - ljduarte@iq.usp.br

This script is made to demonstrate the use of the spectral representation 
of the t and F tests in a 2^2 factorial design, as discussed in the paper: 
Determination of individual metabolic abundance changes owing to environmental 
impacts: t and F-distribution factorial design spectral representations. By 
Gustavo G. Marcheafave, Leonardo J. Duarte, Elis D. Pauli, Ieda S. Scarminio and
Roy E. Bruns.
"""

# Import Libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

"""
 Define functions to calculate F-Distributions
"""
def factorial(a):
    if a > 1:
        value = a * factorial(a-1)
        return value
    if a <= 1:
        return 1   
    
def Beta_function(m, n):
    value = factorial(m-1)*factorial(n-1)/factorial(m+n-1)
    return value

def F_distribution(x, d1, d2):
    value = (1/Beta_function(d1/2, d2/2)) * ((d1/d2)**(d1/2)) * (x**(d1/2 - 1)) * (1 + (d1/d2)* x)**((-d1-d2)/2)
    return value

def Find_interval(d1, d2, alpha):
    x = [i/100 for i in range(1, 10000000)]
    y = [F_distribution(x[i], d1, d2) for i in range(len(x))]
    
    total_area = 0
    for i in range(1, len(x)):
        total_area = total_area +  (y[i-1]+y[i])*(x[i]-x[i-1])*0.5
    
    partial_area = 0
    for i in range(1, len(x)):
        partial_area = partial_area  + (y[i-1]+y[i])*(x[i]-x[i-1])*0.5 
        if partial_area/total_area >= alpha:
            return x[i]

"""
Define function to generate F Spectrum plots.
"""
def F_spectrum(x_values, f_10, f_05, f_01, X, X_e, output_name): 
    f_values_plot = []
    for i in range(X.shape[1]):
        a = sum(np.array(X[:, i])**2)/d1
        b = sum(np.array(X_e[:, i])**2)/d2
        f_values_plot.append(a[0]/b[0])
    
    original_stdout = sys.stdout # Save a reference to the original standard output

    with open(output_name + '.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        for i in range(len(f_values_plot)):
            print('{0:.5f};{1:.5f}'.format(x_values[i],f_values_plot[i]))
        sys.stdout = original_stdout # Reset the standard output to its original value
    count_10=0
    count_05=0
    count_01=0
    
    for value in f_values_plot:
       if value >= f_10:    
            count_10 = count_10 +1
       if value >= f_05:
            count_05 = count_05 +1
       if value >= f_01:
            count_01 = count_01 +1
            
    fig, ax = plt.subplots(1)
    ax.plot(x_values, f_values_plot, color='#8DAA9D', alpha = 1)
    ax.axhline(f_10, ls='--', color='red', label='90%' )
    ax.axhline(f_05, ls='--', color='blue',  label='95%' )
    ax.axhline(f_01, ls='--', color='orange', label='99%' )
    ax.legend()
    ax.set_xlabel('', fontsize=12)
    ax.set_ylabel('F value', fontsize=12)
    ax.set_ylim(0, 80)
    ax.text(0.02, 0.9, 'Significant variables (90%): ' + str(count_10) , size=12, transform=ax.transAxes)
    ax.text(0.02, 0.85, 'Significant variables (95%): ' + str(count_05) , size=12, transform=ax.transAxes)
    ax.text(0.02, 0.80, 'Significant variables (99%): ' + str(count_01) , size=12, transform=ax.transAxes)
    ax.grid(zorder =-500, alpha = 0.5)


    fig.savefig(output_name + '.png', dpi=300)
   
    return None

"""
Define functions to calculate effects and its variances
Xp = Data matrix where each line contains a response variable and each column
contains one spectrum.The first column of XP is the x dimension of the spectrum.
"""

def DOE(Xp, D): #Design of Experiments
    A = np.linalg.inv(np.matmul(D.transpose(), D))
    B = np.matmul(D.transpose(), Xp)
    S = np.matmul(A, B)
    
    mean = np.matrix([[1, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0],
                      [0, 0, 0, 0]])
    
    e1 = np.matrix([[0, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0]])
    
    e2 = np.matrix([[0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 0]])
    
    e12 = np.matrix([[0, 0, 0, 0],
                     [0, 0, 0, 0],
                     [0, 0, 0, 0],
                     [0, 0, 0, 1]])

    X_e1 = np.matmul(np.matmul(D, e1), S)
    X_e2 = np.matmul(np.matmul(D, e2), S)
    X_e12 = np.matmul(np.matmul(D, e12), S)
    X_mean = np.matmul(np.matmul(D, mean), S)

    #Prepare the effect matrices 
    efeitos = np.matmul(D.T, Xp)/D.shape[0]
    efeitos[1:, :] = efeitos[1:, :]*2
    
    X_e = Xp - X_mean - X_e1 - X_e2 - X_e12
    denominator = np.linalg.norm(Xp)**2 - np.linalg.norm(X_mean)**2    
    var_e1 = 100*np.linalg.norm(X_e1)**2/denominator
    var_e2 = 100*np.linalg.norm(X_e2)**2/denominator
    var_e12 =100* np.linalg.norm(X_e12)**2/denominator
    var_e = 100*np.linalg.norm(X_e)**2/denominator
    
    return var_e1, var_e2, var_e12, var_e, X_e1, X_e2, X_e12, X_mean, X_e, efeitos

"""
Function to perform the T-test in all the variables of Xp
"""
def test_T(Xp, N, efeitos):
    S2_list = []
    for i in range(0, 4*N, N):
        m_bloco = Xp[i:i+N,:].sum(axis=0)/(N)
        S2_bloco = ((Xp[i:i+N,:] - m_bloco)**2).sum(axis=0)/(N-1)
        S2_list.append(S2_bloco)
    Ep = (sum(S2_list)/(4*N))**0.5
    T = efeitos[1:, :]
    for i in range(len(Ep)):
        T[:,i] = T[:,i]/Ep[i]

    return Ep, T    

"""
Define function to generate histogram
"""
def Plot_hist(fe_values, p_value, output_name):  
    fig, ax = plt.subplots(1)
    n, bins, patches = ax.hist(fe_values, 50, density=False, facecolor='b', edgecolor='k', linewidth=1, alpha=0.75)
    ax.axvline(fe_values[0], c='r')
    ax.set_ylabel('Counts', size=12)
    ax.set_xlabel('% Variance', size=12)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.text(0.7, 0.9, 'p-value: ' + str(p_value) , size=12, transform=ax.transAxes)
    fig.savefig(output_name + '.png', dpi=300)
    plt.clf()
    
    return None

"""
Generate the T-spectrum figure.
"""
def plot_T(T_values, x_values, N, output_name):
    #N=number of replicates
    lim = 0
    if N == 2:
        lim = 2.776
    elif N== 3:
        lim = 2.306
    elif N == 4:
        lim = 2.179
    elif N == 5:
        lim = 2.120
    titles = ['Factor 1', 'Factor 2', 'Interaction 12']
        
    for i in range(T_values[1].shape[0]):
        original_stdout = sys.stdout
        with open(output_name + '_T_' + titles[i] + '.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            for j in range(len(T_values[1][i, :])):
                print('{0:.5f};{1:.5f}'.format(x_values[j],T_values[1][i, :][j]))
            sys.stdout = original_stdout # Reset the standard output to its original value
        fig, ax = plt.subplots(1)
        ax.set_title(titles[i])
        ax.plot(x_values, T_values[1][i, :], color="blue", alpha=0.8, zorder=500)
        ax.axhline(lim, c='r', lw=1, ls="--", zorder=-500)
        ax.axhline(-1*lim, c='r', lw=1, ls='--', zorder=-500)
        ax.axhline(0, c='k', lw=1, ls='--', zorder=-500)
        ax.grid(zorder=-500)
        ax.set_xlim(min(x_values), max(x_values))
        ax.set_ylim(-max(abs(T_values[1][i, :]))-0.2, max(abs(T_values[1][i, :]))+0.2)
        ax.set_xlabel("")
        ax.set_ylabel("T Value")
        ax.locator_params(axis='both', nbins=10)
        fig.savefig(output_name + '_T_' + titles[i] +'.png', dpi=300)
    return None


if __name__ == "__main__":
     #ler os dados de input:
     input_name = input("Enter the data matrix in the .csv format:")
     DOE_matrix_name = input("Insert the DOE matrix in the .csv format: ")
     output_name = input_name[0:-4] + "_output"
     
     """
     Calculate the DOE effect values
     """
     #### read the data matrix
     data = pd.read_csv(input_name, header=None, sep=';')
     x_values = data[0]
     Xp = np.array(data.loc[:, 1:].transpose()) # response
     #### read the factorial matrix
     D = pd.read_csv(DOE_matrix_name, header=None, sep=';')
     D = np.array(D)
     ### calculate the effects:
     var_e1, var_e2, var_e12, var_e, X_e1, X_e2, X_e12, X_mean, X_e, efeitos = DOE(Xp, D)
        
     ### report
     print('########## Effects variance: ###########')
     print('% Variance of E1: ',var_e1)
     print('% Variance of E2: ',var_e2)
     print('% Variance of E12: ',var_e12)
     print('% Variance of Residual: ',var_e)
     print('% Total: ',var_e1 + var_e2 + var_e12 + var_e)
    
     """
     Perform the t test in each variable and generate figures
     """
     N = int(D.shape[0]/4)
     T_values = test_T(Xp, N, efeitos)
     plot_T(T_values, x_values, N, output_name)
     
     """
     Perform the F test in each variable and generate figures
     """
     d2 = D.shape[0]-4
     d1 = 1
     f_10 = Find_interval(d1, d2, 0.90)
     f_05 = Find_interval(d1, d2, 0.95)
     f_01 = Find_interval(d1, d2, 0.99)
     F_spectrum(x_values, f_10, f_05, f_01, X_e1, X_e, output_name + '_F_dist_e1')
     F_spectrum(x_values, f_10, f_05, f_01, X_e2, X_e, output_name + '_F_dist_e2')
     F_spectrum(x_values, f_10, f_05, f_01, X_e12, X_e, output_name + '_F_dist_e12')
     
     """
     Print the results
     """
     original_stdout = sys.stdout # Save a reference to the original standard output

     with open(output_name + '_e1_DOE_.txt', 'w') as f:
         sys.stdout = f # Change the standard output to the file we created.
         print_e1 = np.array(X_e1[0])[0]
         for i in range(len(x_values)):
             print("{0:2.5f};{1:2.5f}".format(x_values[i], print_e1[i]))
     sys.stdout = original_stdout # Reset the standard output to its original value
     
     with open(output_name + '_e2_DOE_.txt', 'w') as f:
         sys.stdout = f # Change the standard output to the file we created.
         print_e2 = np.array(X_e2[0])[0]
         for i in range(len(x_values)):
             print("{0:2.5f};{1:2.5f}".format(x_values[i], print_e2[i]))
     sys.stdout = original_stdout # Reset the standard output to its original value
     
     with open(output_name + '_e12_DOE_.txt', 'w') as f:
         sys.stdout = f # Change the standard output to the file we created.
         print_e12 = np.array(X_e12[0])[0]
         for i in range(len(x_values)):
             print("{0:2.5f};{1:2.5f}".format(x_values[i], print_e12[i]))
     sys.stdout = original_stdout # Reset the standard output to its original value
     
     


