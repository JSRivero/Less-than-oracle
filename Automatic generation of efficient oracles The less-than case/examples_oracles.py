# Imports

import functions_oracles as func


# Some examples:

if __name__=='__main__':

    nqubits = 6

    # Less than oracle
    number_less_than = 42
    less_than_oracle = func.oracle_less_than(number=number_less_than, nqubits=nqubits)

    # Greater than oracle
    number_larger_than = 23
    greater_than_oracle = func.oracle_greater_than(number=number_larger_than, nqubits=nqubits)

    # Interval oracle
    lower_boundary = 25
    upper_boundary = 42
    interval_oracle = func.oracle_interval(lower_boundary=lower_boundary, upper_boundary=upper_boundary, nqubits=nqubits)

    #########################
    #### Unitary oracles ####
    #########################

    # Less than oracle
    number_less_than = 42
    unitary_less_than_oracle = func.diagonal_oracle_less_than(number=number_less_than, nqubits=nqubits, name='Diagonal Oracle')
