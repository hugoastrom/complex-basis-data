#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import math
from matplotlib.ticker import StrMethodFormatter


bases = ['cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ', 'cc-pV5Z','aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'aug-cc-pV5Z', 'HGBSP1-5', 'HGBSP1-7', 'HGBSP1-9', 'HGBSP2-5', 'HGBSP2-7', 'HGBSP2-9', 'HGBSP3-5', 'HGBSP3-7', 'HGBSP3-9', 'AHGBSP1-5', 'AHGBSP1-7', 'AHGBSP1-9', 'AHGBSP2-5', 'AHGBSP2-7', 'AHGBSP2-9', 'AHGBSP3-5', 'AHGBSP3-7', 'AHGBSP3-9']#, '6-311++G(3df,3pd)', 'def2-TZVP']

#atoms = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar']
atoms = ["He"]

# Recreate figures?
do_figures = False

def load_gto_energy(basis, atoms):
    '''Load GTO energies from Erkale log file'''
    data = {}
    for at in atoms:
        states = {}
        for state in range(n_states(at)):
            mag_fields = {}
            for field in np.arange(0.00, 0.62, 0.02):
                field = '{:.2f}'.format(field)
                try:
                    f = open(f'output/{at}/{field}/{basis}_{state}.log')
                except:
                    mag_fields[field] = None
                for line in f:
                    if 'Total energy:' in line:
                        line_split = line.split()
                        mag_fields[field]=float(line_split[-1])
                        break
            states[state]=mag_fields
        data[at]=states

    return data

#def load_complex_energy(basis, atoms):
#    '''Load zero field energy from QChem output file with complex orbitals'''
#    data = {}
#    for at in atoms:
#        try:
#            f = open(f'{path}')
#        except:
#            data[at] = None
#        for line in f:
#            if 'Total energy:' #input correct line later
#                line_split = line.split()
#                data[at] = float(line_split[-1])
#                break

#    return data


def n_states(at):
    n=0
    while os.path.exists(f'configs/{at}_{n}.occs'):
        n+=1
    return n

def load_fem_energy(atoms):
    '''Load FEM energies from .log files'''
    data = {}
    for at in atoms:
        states = {}
        for state in range(n_states(at)):
            mag_fields = {}
            for field in np.arange(0.00, 0.70, 0.10):
                field = '{:.2f}'.format(field)
                lmax = 21
                while True:
                    if os.path.exists(f'output/{at}/{field}/fem_{lmax}_{state}.log'):
                        break
                    lmax = lmax - 1
                f = open(f'output/{at}/{field}/fem_{lmax}_{state}.log')
                for line in f:
                    if 'Total                 energy:' in line:
                        line_split=line.split()
                        mag_fields[field]=float(line_split[-1])
                        break
            states[state]=mag_fields
        data[at]=states

    return data

def difference_is_positive(atom,state):
    '''check if bste is positive'''
    bste_positive = True
    for field in all_results['FEM'][atom][state]:
        dE = all_results['AHGBSP3-9'][atom][state][field] - all_results['FEM'][atom][state][field]
        if dE < -1e-6:
            bste_positive = False
    return bste_positive


def compute_difference_vector(basis,atom,state):
    try:
        return [ all_results[basis][atom][state][field] - all_results['FEM'][atom][state][field] for field in all_results['FEM'][atom][state] ]
    except:
        return None


def compute_avg_bste(basis,atom,state):
    error_vector = compute_difference_vector(basis,atom,state)
    if error_vector is not None:
        return np.mean(np.absolute(error_vector))
    else:
        return None


def make_labels(at):
    '''make state labels'''
    occ_list = []
    for i in range(n_states(at)):
        occ = '$'

        # n=0...3 means sigma...phi orbitals
        for n in range(4):
            occfile=open(f'configs/{at}_{i}.occs')
            for line in occfile:
                line_split = line.split()
                if int(line_split[-1]) == n or int(line_split[-1]) == -n:
                    occ+=mapp[line_split[-1]]
                    occ+=f'^{{{line_split[0]},{line_split[1]}}}'
        occ+='$'
        occ_list.append(occ)
    return occ_list



def table_one_atom(atom,subset,fname,label,caption,sideways):
    '''Generate LaTeX table with BSTE for atom'''
    # open .tex file
    texfile=open('tables/'+fname,'w')
    texfile.write('\\begin{{{}*}}\n'.format('sidewaystable' if sideways else 'table'))
    texfile.write('\centering\n')
    texfile.write('\\tiny\n')
    texfile.write('\\begin{{tabular}}{{l{}}}\n'.format('c'*n_states(atom)))
    texfile.write('\hline\n')

    # header
    for state in make_labels(atom):
        texfile.write(f' & {state}')
    texfile.write('\\\\ \n')
    texfile.write('\hline \hline\n')

    # write data
    for basis in subset:
        texfile.write(basis)
        for state in range(n_states(atom)):
            if compute_avg_bste(basis,atom,state) is not None:
                color = 'blue' if difference_is_positive(atom,state) else 'red'
                texfile.write(' & \color{{{}}}{{$ {:.2f} $}}'.format(color, 1000*compute_avg_bste(basis,atom,state)))
            else:
                texfile.write(' & ')
        texfile.write('\\\\ \n')
    texfile.write('\hline\n')
    texfile.write('\end{tabular}\n')
    texfile.write(f'\caption{{{caption}}}\n')
    texfile.write(f'\label{{tab:{label}}}\n')
    texfile.write('\\end{{{}*}}\n'.format('sidewaystable' if sideways else 'table'))
    texfile.close()

def table_one_basis(atoms,basis,fname,label,caption,sideways):
    '''Generate LaTeX table with BSTE for one basis'''
    # open .tex file
    texfile=open('tables/'+fname,'w')
    texfile.write('\\begin{{{}*}}\n'.format('sidewaystable' if sideways else 'table'))
    texfile.write('\centering\n')
    texfile.write('\\tiny\n')
    texfile.write('\\begin{{tabular}}{{l{}}}\n'.format('c'*9))
    texfile.write('\hline\n')

    # header
    for i in range(9):
        texfile.write(f' & state {i} ')
    texfile.write('\\\\ \n')
    texfile.write('\hline \hline\n')

    # write data
    for at in atoms:
        texfile.write(at)
        for state in make_labels(at):
            texfile.write(f' & {state}')
        texfile.write(' &'*(9-n_states(at)))
        texfile.write('\\\\ \n')
        for state in range(n_states(at)):
            if compute_avg_bste(basis,at,state) is not None:
                color = 'blue' if difference_is_positive(at,state) else 'red'
                texfile.write(' & \color{{{}}}{{$ {:.2f} $}}'.format(color, 1000*compute_avg_bste(basis,at,state)))
            else:
                texfile.write(' & ')
        texfile.write(' &'*(9-n_states(at)))
        texfile.write('\\\\ \n')
    texfile.write('\hline\n')
    texfile.write('\end{tabular}\n')
    texfile.write(f'\caption{{{caption}}}\n')
    texfile.write(f'\label{{tab:{label}}}\n')
    texfile.write('\\end{{{}*}}\n'.format('sidewaystable' if sideways else 'table'))
    texfile.close()

def table_total_energy(basis, atom, fname, label, caption):
    '''Generate LaTeX table with total energies for all field strengths'''
    # open .tex file
    texfile=open('tables/'+fname,'w')
    texfile.write('\\begin{table*}\n')
    texfile.write('\centering\n')
    texfile.write('\\tiny\n')

    # Field strengths to print out
    fields = [field for field in all_results['FEM'][atom][0] if field in all_results[basis][atom][0]]

    texfile.write('\\begin{{tabular}}{{l{}}}\n'.format('c'*len(fields)))
    texfile.write('\hline\n')

    # header
    for field in fields:
        texfile.write(f' & $ {field} B_0 $')
    texfile.write('\\\\ \n')
    texfile.write('\hline \hline\n')

    # write data
    for state in range(n_states(atom)):
        texfile.write(make_labels(atom)[state])
        for field in fields:
            if field in all_results[basis][atom][state] and all_results[basis][atom][state][field] is not None:
                # Energies are converged to 1 uEh so report to that accuracy
                texfile.write(' & $ {:.6f} $'.format(all_results[basis][atom][state][field]))
            else:
                texfile.write(' & ')
        texfile.write('\\\\ \n')

    texfile.write('\hline\n')
    texfile.write('\end{tabular}\n')
    texfile.write(f'\caption{{{caption}}}\n')
    texfile.write(f'\label{{tab:{label}}}\n'.replace(',',''))
    texfile.write('\end{table*}\n')
    texfile.close()

def table_mean_difference(bases, atom, fname, label, caption):
    '''Generate LaTeX table with mean absolute energy differences between FEM and GTO for all states'''
    # open .tex file
    texfile=open('tables/'+fname,'w')
    texfile.write('\\begin{table}\n')
    texfile.write('\centering\n')
    texfile.write('\\small\n')
    texfile.write('\\begin{{tabular}}{{ll{}}}\n'.format('c'*len(bases)))
    texfile.write('\hline\n')

    # header
    texfile.write(' & state')
    for basis in bases:
        texfile.write(f' & {basis}')
    texfile.write('\\\\ \n')
    texfile.write('\hline \hline\n')

    # write data
    for state in range(n_states(atom)):
        texfile.write(f'{state} & {make_labels(atom)[state]}')
        for basis in bases:
            texfile.write(' & ')
            diff_vector = compute_difference_vector(basis,atom,state)
            if diff_vector is not None:
                diff = np.mean(np.absolute(diff_vector))
                color = 'blue' if difference_is_positive(atom,state) else 'red'
                texfile.write('\color{{{}}}{{ $ {:.3f} $ }}'.format(color, 1000*diff))
        texfile.write('\\\\ \n')

    texfile.write('\hline\n')
    texfile.write('\end{tabular}\n')
    texfile.write(f'\caption{{{caption}}}\n')
    texfile.write(f'\label{{tab:{label}}}\n')
    texfile.write('\end{table}\n')
    texfile.close()

def find_min_max(at):
    '''finds min and max FEM energies of atom and returns y-axis limits for figures'''
    Emin=0.0
    Emax=-1000.0

    # find min and max FEM energy
    for state in range(n_states(at)):
        for field in all_results['FEM'][at][state]:
            if Emin>all_results['FEM'][at][state][field]:
                Emin=all_results['FEM'][at][state][field]
            if Emax<all_results['FEM'][at][state][field]:
                Emax=all_results['FEM'][at][state][field]

    # find largest difference between FEM and aug-cc-pVTZ
    diff = 0.0
    for state in range(n_states(at)):
        for field in all_results['FEM'][at][state]:
            if diff < all_results['aug-cc-pVTZ'][at][state][field] - all_results['FEM'][at][state][field]:
                diff = all_results['aug-cc-pVTZ'][at][state][field] - all_results['FEM'][at][state][field]
    
    Emax += diff/2 + 0.25
    return Emin, Emax

def plot_energies(at,basis):
    '''plots energy as a function of magnetic field strength for atom in basis'''

    for state in range(n_states(at)):
        if  all_results[basis][at][state] is not None:

            # FEM B fields and energies
            b_fields= [ float(field) for field in all_results['FEM'][at][state] ]
            fem_energies = [ all_results['FEM'][at][state][field] for field in all_results['FEM'][at][state] ]
            plt.plot(b_fields, fem_energies, color = colors[state], linewidth=0.0, marker='s', label=f'{state}')

            # GTO B fields and energies
            b_fields= [ float(field) for field in all_results[basis][at][state] ]
            gto_energies = [ all_results[basis][at][state][field] for field in all_results[basis][at][state] ]
            plt.plot(b_fields, gto_energies, color = colors[state])

def plot_convergence(at):
    '''plots convergence to CBS limit of numerical basis'''

    maxl=0
    labels = make_labels(at)

    for state in range(n_states(at)):
        lvals = []
        energy = []

        # assign start l value
        f = open(f'configs/{at}_{state}.occs')
        mvals = [ int(line.split()[-1]) for line in f ]
        l = max(np.absolute(mvals))

        # append l values and energy of that l value
        while os.path.exists(f'output/{at}/0.60/fem_{l}_{state}.log'):
            lvals.append(l)
            f = open(f'output/{at}/0.60/fem_{l}_{state}.log')
            for line in f:
                if 'Total                 energy:' in line:
                    line_split=line.split()
                    energy.append(float(line_split[-1]))
                    break
            l+=2

        # find maximum l value and setup plot
        if int(lvals[-1])>maxl:
            maxl = int(lvals[-1])
        diff = [ energy[i] - energy[i+1] for i in range(len(lvals)-1) ]
        del lvals[0]
        plt.semilogy(lvals, diff, color=colors[state], linewidth=0.0, marker='s', label=labels[state])

    # convergence limit 1 uEh
    plt.plot(np.arange(maxl+1), 1e-6*np.ones(maxl+1), 'k--')



def print_fig(at, basis, fname, xlabel, ylabel, energy_plot):
    '''generate figures'''
    # figure parameters
    fig=plt.figure()
    fig.set_figheight(10)
    fig.set_figwidth(13.5)
    figsize=(6,6)
    fontsize=23
    # one decimal in y axis values
    plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    if energy_plot:
        Emin = find_min_max(at)[0]
        Emax = find_min_max(at)[1]
        plot_energies(at, basis)
        Emin=(math.floor(Emin*10))/10
        Emax=(math.ceil(Emax*10))/10
        plt.ylim(Emin, Emax)
        plt.legend(loc=9, ncol=3, fontsize=fontsize)
    else:
        plot_convergence(at)
        plt.legend(loc='upper right', fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    plt.savefig(f'figures/{fname}.png')
    plt.close()


def print_toc(atom_subset, basis_subset, fname, ylabel):
    '''Print TOC violin plot'''
    fontsize=23
    errors = []
    xlabel = []
    for atom in atom_subset:
        for basis in basis_subset:
            ers_basis = [ compute_avg_bste(basis,atom,state) for state in range(n_states(atom)) ]
            errors.append(ers_basis)
            xlabel.append(f'{atom} {basis}')

    fig, ax = plt.subplots(figsize=(13, 7))
    plots = ax.violinplot(errors, vert=True, showmedians=True, showextrema=True, widths=1)

    # set colors of violin patches
    for pc, color in zip(plots['bodies'], colors):
        pc.set_facecolor(color)

    # set color of median line
    for param in ['cmedians', 'cmins', 'cmaxes', 'cbars']:
        plots[param].set_colors(color)

    # set labels
    num_plots = [ i+1 for i in range(len(errors)) ]
    ax.set_xticks(num_plots, labels=xlabel, fontsize=fontsize, rotation=90)
    ax.set_ylabel(ylabel, fontsize=fontsize)

    plt.axhline(color='grey', linestyle='-')
    plt.tight_layout()
    plt.savefig(f'{fname}.png')

#def zero_field_table(realdata, complexdata, basis, fname, label, caption):
#    '''Print table of ground state energies at zero field'''

    # open latex file
#    texfile=open(f'tables/{fname}','w')
#    texfile.write('\\begin{table}\n')
#    texfile.write('\centering\n')
#    texfile.write('\\small\n')
#    texfile.write('\\begin{{tabular}}{{clll}}\n')
#    texfile.write('\hline\n')

    # header
#    texfile.write(' & Real & Complex & Numerical \\\\ \n')
#    texfile.write('\hline \hline\n')

#    for at in atoms:
#        texfile.write(f'{at}')
#        texfile.write(f'{realdata} & {complexdata}')

mapp={
    '-3': '\phi_+',
    '-2': '\delta_+',
    '-1': '\pi_+',
    '0': '\sigma',
    '1': '\pi_-',
    '2': '\delta_-',
    '3': '\phi_-'
}

colors = ['b','tab:orange','g','r','c','tab:pink','y','k','tab:brown']


# load all results
all_results = {'FEM' : load_fem_energy(atoms)}
for basis in bases:
    all_results[basis] = load_gto_energy(basis, atoms)


# write results to csv files
for basis in all_results:
    for at in all_results[basis]:
        csvfile = open(f'csvfiles/{at}_{basis}.csv','w')
        for field in all_results[basis][at][0]:
            csvfile.write(field)
            for state in all_results[basis][at]:
                try:
                    csvfile.write(f',{all_results[basis][at][state][field]}')
                except:
                    csvfile.write(',None')
            csvfile.write('\n')
        csvfile.close()


# produce tables
for at in atoms:
    table_one_atom(at, bases, f'{at}_gto_tab.tex', f'{at}-gto', f'Mean absolute energy differences $\Delta E^\\text{{GTO}}$ in m$E_h$ for a variety of GTO basis sets for the {at} atom. States with all positive differences $\Delta E$ are shown in blue, and states that exhibit negative $\Delta E$ at one or more field strength(s) in red.', sideways=(n_states(at) > 6))

for at in atoms:
    table_total_energy('FEM', at, f'{at}_fem_tab.tex', f'{at}-fem', f'Complete basis set limit total energies in $E_h$ for the {at} atom at all field strengths.')
    for basis in bases:
        table_total_energy(basis, at, f'{at}_{basis}_tab.tex', f'{at}-{basis}', f'Total energies in $E_h$ for the {at} atom in the {basis} basis set in fully uncontracted form, employing the real-orbital approximation.')

subset = ['aug-cc-pVTZ', 'AHGBSP3-9']

for basis in subset:
    table_one_basis(atoms,basis,f'{basis}_tab.tex',f'{basis}',f'Energy differences $\Delta E$ in m$E_h$ for the fully uncontracted {basis} basis set. States with all positive $\Delta E$ shown in blue, and states that exhibit negative $\Delta E$ at one or more field strength shown in red.', sideways = True)

for at in atoms:
    table_mean_difference(subset, at, f'{at}-mean-diff.tex', f'{at}-mean-differ', f'MAEDs between GTO and FEM energies in m$E_h$ for {at} in the fully uncontracted {" and ".join(subset)} basis sets.')


# LaTeX input file for the manuscript figures
for at in atoms:
    figfile = open('figures/manuscript-{}.tex'.format(at), 'w')
    figfile.write('\\begin{figure*}\n')
    figfile.write('\\begin{center}\n')
    for basis in subset:
        figfile.write(f'\includegraphics[width=0.45\linewidth]{{figures/{at}-{basis}.png}}\n')
    figfile.write(f'\caption{{Total energy of the {at} atom as a function of the magnetic field strength $B$ in the {subset[0]} (left) and {subset[1]} (right) basis sets.}}\n')
    figfile.write(f'\label{{fig:{at}}}\n')
    figfile.write('\end{center}\n')
    figfile.write('\end{figure*}\n')
    figfile.write('\n')
    figfile.close()

# plots of convergence of numerical basis to CBS limit as function of truncation parameter L
figures = open('figures/convergence-figures.tex', 'w')
for at in atoms:
    figures.write('\\begin{figure}\n')
    figures.write(f'\includegraphics[width=\linewidth]{{figures/convergence_{at}.png}}\n')
    figures.write(f'\caption{{Convergence of the total energy of the considered states of the {at} atom as a function of the maximum angular momentum included in the fully numerical basis set.}}\n')
    figures.write(f'\label{{fig:{at}}}\n')
    figures.write('\end{figure}\n')
    if do_figures:
        print_fig(at, 'fem', f'convergence_{at}', 'L', '$\Delta E$', False)
figures.close()

# plots of GTO BSTEs as function of magnetic field B
figures = open('figures/figures.tex', 'w')
for at in atoms:
    for basis in bases:
        # energy plots
        figures.write('\\begin{figure}[H]\n')
        figures.write(f'\includegraphics[width=\linewidth]{{figures/{at}-{basis}.png}}\n')
        figures.write(f'\caption{{Total energies of all considered states of the {at} atom in the {basis} basis set in fully uncontracted form (solid lines). The FEM values are shown by the squares of the same color.}}\n')
        figures.write(f'\label{{fig:{at}-{basis}}}\n'.replace(',',''))
        figures.write('\end{figure}\n')
        if do_figures:
            print_fig(at, basis, f'{at}-{basis}', 'B/$B_0$', 'E/$E_h$', True)
    figures.write('\clearpage\n')


figures.close()

# generate SI text
SItext = open('manuscript/SItext.tex', 'w')

SItext.write('\section{Convergence of Total Energies in FEM to the CBS Limit}')
SItext.write('We begin by showing that the complete basis set (CBS) limit is achieved with the employed finite element method (FEM) by plotting the convergence of the total energy at the field strength $B=0.60$ as a function of the angular truncation parameter $l_\\text{max}$.\n')
SItext.write('The energy difference $\\Delta E = E(l_\\text{max}) - E(l_\\text{max}-2)$ is shown\n')
for iat, at in enumerate(atoms):
    if iat>0:
        SItext.write(', ')
    if iat == len(atoms)-1:
        SItext.write('and ')
    SItext.write(f'in \\cref{{fig:{at}}} for {at}')
SItext.write('.\n\n')

SItext.write('The resulting CBS total energies of all atoms, determined with FEM, are given in\n')
for iat, at in enumerate(atoms):
    if iat>0:
        SItext.write(', ')
    if iat == len(atoms)-1:
        SItext.write('and ')
    SItext.write(f'in \\cref{{tab:{at}-fem}} for {at}')
SItext.write('.\n\n')

SItext.write('\section{Difference of GTO Energies to FEM Values}')
SItext.write('\subsection{Mean Absolute Differences}')
SItext.write('Given these total energies, the differences $\Delta E^\\text{GTO} = E^\\text{GTO} - E^\\text{CBS}$ at each field strength for each state are obtained.\n')
SItext.write('First, the average absolute energy differences for the studied Gaussian basis sets in the fully uncontracted form are given in\n')
for iat, at in enumerate(atoms):
    if iat>0:
        SItext.write(', ')
    if iat == len(atoms)-1:
        SItext.write('and ')
    SItext.write(f'in \\cref{{tab:{at}-gto}} for {at}')
SItext.write('.\n\n')

SItext.write('Missing entries in the tables indicate either that the basis set for the given element does not exist on the Basis Set Exchange, e.g. aug-cc-pV5Z for Li in \\cref{tab:Li-gto}, or that the basis set is too small to describe the state in question, e.g. the state of the C atom with the occupied $\\varphi$ orbital in \\cref{tab:C-gto} which requires at least $f$ functions in the atomic basis.\n\n')

SItext.write('\subsection{Plots of FEM and GTO Total Energies}')
SItext.write('Next, we include plots of the differences for all the studied states of all the studied atoms in all the studied basis sets as a function of the field strength $B$.\n\n')
for iat, at in enumerate(atoms):
    SItext.write(f'The results for the {at} atom are given in\n')
    for ibasis, basis in enumerate(bases):
        if ibasis>0:
            SItext.write(', ')
        if ibasis == len(bases)-1:
            SItext.write('and ')
        SItext.write(f'in \\cref{{fig:{at}-{basis}}} for {basis}'.replace(',',''))
    SItext.write(', all basis sets being employed in the fully uncontracted form.\n\n')

SItext.write('\subsection{Tables of GTO Total Energies}')
SItext.write('Finally, we report the state specific total energies for all the studied states of all the studied atoms in all the studied Gaussian basis sets as a function of the field strength $B$, employing the real-orbital approximation.\n\n')
for at in atoms:
    SItext.write(f'For {at}, the results are given in\n')
    for ibasis, basis in enumerate(bases):
        if ibasis>0:
            SItext.write(', ')
        if ibasis == len(bases)-1:
            SItext.write('and ')
        SItext.write(f'in \\cref{{tab:{at}-{basis}}} for {basis}'.replace(',',''))
    SItext.write(', all basis sets being employed in the fully uncontracted form.\n\n')


# Input the corresponding tables and figures here
SItext.write('\clearpage\n')
SItext.write('\input{figures/convergence-figures.tex}\n')
SItext.write('\clearpage\n')
for at in atoms:
    SItext.write(f'\\input{{tables/{at}_fem_tab.tex}}\n')
SItext.write('\n')

SItext.write('\clearpage\n')
for at in atoms:
    SItext.write(f'\\input{{tables/{at}_gto_tab.tex}}\n')
SItext.write('\clearpage\n')
SItext.write('\input{figures/figures.tex}\n')
SItext.write('\n')

SItext.write('\clearpage\n')
for at in atoms:
    for basis in bases:
        SItext.write(f'\\input{{tables/{at}_{basis}_tab}}\n')
    SItext.write('\clearpage')
SItext.write('\n')
SItext.close()

# print TOC figure
if do_figures:
    print_toc(['He', 'F', 'Si', 'Ar'], ['aug-cc-pVTZ', 'AHGBSP3-9'], 'figures/toc', 'Difference from FEM [$E_h$]')
