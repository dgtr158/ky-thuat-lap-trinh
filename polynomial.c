#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Polynomial abstract data type
typedef struct polynomial {
    int degree;
    double* coefficient; 
} Polynomial;


/* ------------ GLOBAL VARIABLES ------------ */
// The list to hold all the intervals containing roots of the equation
double* rootIntervals;

// The file for storing the output of the program
FILE *output;

/* ------------ FUNCTION'S PROTOTYPE ------------ */
// Function to display menu
void displayMenu();

// Function to input the polynomial
Polynomial initPolynomial();

// Function to display interval
void displayInterval(int size);

// Function to get the value of polynomial at point x
double getValue(Polynomial polynomial, double x);

// Function to get the first derivative function of a polynomial 
Polynomial getDerivativeFunction(Polynomial polynomial);

// Function to find upper bound of range of real roots
double findUpperBound(Polynomial polynomial);

// Function to find lower bound of range of real roots
double findLowerBound(Polynomial polynomial);

// Function to shrink an interval containing root
double* shrinkIntervalContainingRoot(Polynomial polynomial, double* interval);

// Function to define error function 1
double firstErrorFunction(Polynomial polynomial, double* interval, double x);

// Function to define error function 2
double secondErrorFunction(Polynomial polynomial, double *interval, double cur, double prev);

// Function to find a root with constant number of iteration
double findRootWithConstantIteration(Polynomial polynomial, double *interval);

// Function to find a root with given error e and use the first error function
double findRootWithGivenErrorUseFirstErrorFunc(Polynomial polynomial, double *interval);

// Function to find a root with given error e and use the second error function
double findRootWithGivenErrorUseSecondErrorFunc(Polynomial polynomial, double *interval);

// Function to find a root with given error e satisfying |Xn - Xn-1| <= e
double findRootWithGivenError(Polynomial polynomial, double *interval);

// Function to check if an interval contains a real root
int hasRealRoot(Polynomial polynomial, double start, double end);

// Gradient decent function
double gradientDescent(Polynomial polynomial, double guess, double upperBound);

// find all interval containing roots
int findAllIntervalContainingRoots(Polynomial polynomial);

// Function to display shrinked interval
void displayShrinkedInterval(double* interval, int size);

/* ------------ MAIN FUNCTION ------------ */
int main() {

    // Create and open output file in the current directory
    output = fopen("output.txt", "w");
    if (!output) {
        printf("Cannot open this file\n");
        return 1;
    }
    else {
        printf("File is opened\n");
    }

    // Polynomial initialization
    Polynomial polynomial;

    // Main flow
    int option;
    int quit = 0;
    int size; // size of interval containing all roots
    double* shrinkedInterval = (double*)malloc(sizeof(double) * 2); 
    double root1, root2, root3, root4; // for storing the roots found

    while (!quit) {
        displayMenu();
        
        // prompt user for option
        scanf("%d", &option);

        switch (option) {
            case 1:
                // system("clear");
                polynomial = initPolynomial();
                break;
            case 2:
                system("clear");
                size = findAllIntervalContainingRoots(polynomial);
                displayInterval(size);
                quit = (size < 1) ? 1 : 0;
                // test
                printf("size = %d\n", size);
                break;
            case 3:
                system("clear");
                shrinkedInterval[0] = rootIntervals[0];
                shrinkedInterval[1] = rootIntervals[1]; 
                shrinkedInterval = shrinkIntervalContainingRoot(polynomial, shrinkedInterval);
                displayShrinkedInterval(shrinkedInterval, 2);
                break;
            case 4:
                system("clear");
                root1 = findRootWithConstantIteration(polynomial, shrinkedInterval);
                break;   
            case 5:
                system("clear");
                root2 = findRootWithGivenErrorUseFirstErrorFunc(polynomial, shrinkedInterval);
                root3 = findRootWithGivenErrorUseSecondErrorFunc(polynomial, shrinkedInterval);       
                break;
            case 6:
                system("clear");
                root4 = findRootWithGivenError(polynomial, shrinkedInterval);
                break;
            case 0:
                quit = 1;
                break;
        }
    }

    // free
    free(rootIntervals);
    free(shrinkedInterval);

    // close file
    fclose(output);
    return 0;
}

/* ------------ FUNCTION'S DETAILS ------------ */
// Function to display menu
void displayMenu()
{
    printf("\t\t\t  =============================MENU================================\n");
    printf("\t\t\t |1. Enter the polynomial\n");
    printf("\t\t\t |2. Find all intervals containing root of the polynomial\n");
    printf("\t\t\t |3. Find interval containing root [a, b] of the polynomial such that |a - b| <= 0.5\n");
    printf("\t\t\t |4. Find the root of the polynomial given a particular number of iteration in the interval(a, b) and the two error functions\n");
    printf("\t\t\t |5. Find the root of the polynomial given a particular error in the interval(a, b) and the two error functions\n");
    printf("\t\t\t |6. Find the root Xn of the polynomial in the interval (a, b) such that |Xn - Xn-1| <= e (e is entered by user)\n");
    printf("\t\t\t |0. Exit the program\n");
    printf("\t\t\t Your option: ");
}

// Function for inputing the polynomial
Polynomial initPolynomial() {
    // write to file
    fprintf(output, " ===== INPUT POLYNOMIAL ====== \n");

    // declaration
    Polynomial polynomial;
    int i; // loop variable

    // prompt user for degree of polynomial
    do {
        printf("Enter the degree of the polynomial: ");
        scanf("%d", &polynomial.degree);
        if (polynomial.degree < 0) {
            printf("The degree of a polynomial must be greater or equal to 0!\n");
        }
    } while (polynomial.degree < 0);

    // write degree to the file
    fprintf(output, "Degree: %d\n", polynomial.degree);

    // allocate memory for storing polynomial's coefficients
    polynomial.coefficient = (double* ) malloc((polynomial.degree + 1) * sizeof(double));

    // prompt user for coefficients of polynomial
    printf("Enter the coefficients of the polynomial: ");
    fprintf(output, "Coefficients: ");
    for (i = 0; i <= polynomial.degree; i++) {
        scanf("%lf", &polynomial.coefficient[i]);
        fprintf(output, "%7.4lf", polynomial.coefficient[i]); // write the coefficient to output file
    }
    fprintf(output, "\n\n");

    return polynomial;
}

// Function to get the value of polynomial at point x
double getValue(Polynomial polynomial, double x) {
    // declaration
    double value = 0;
    int i; // loop variable

    // get value at point x
    for (i = 0; i <= polynomial.degree; i++) {
        value += polynomial.coefficient[i] * pow(x, (double)(polynomial.degree-i));
    }

    // return 
    return value;
}

// Function to get the first derivative function of a polynomial
Polynomial getDerivativeFunction(Polynomial polynomial) {
    // Variable declaration
    Polynomial df;
    int i; // loop variable

    // get the derivative function
    df.degree = polynomial.degree - 1; // degree of derivative function
    df.coefficient = (double*) malloc (df.degree * sizeof(double));
    for (i = 0; i <= df.degree; i++) {
        df.coefficient[i] = (polynomial.degree-i) * polynomial.coefficient[i];
    }

    // return 
    return df;
}

// Function to shrink an interval containing root
double *shrinkIntervalContainingRoot(Polynomial polynomial, double *interval) {
    // declaration
    double tolerance = 0.5; // restrict the length of interval
    double start = interval[0]; // the begining of interval
    double end = interval[1]; // the end of the interval
    double mid; // midpoint of start and end
    double *shrinkedInterval = (double*) malloc (2 * sizeof(double));

        // shrink the interval until |start - end| <= tolerant
        while (fabs(start - end) > tolerance)
    {
        mid = (start + end) / 2;
        if (getValue(polynomial, mid) * getValue(polynomial, start) > 0) {
            start = mid;
        }
        else {
            end = mid;
        }
    }

    // store interval
    shrinkedInterval[0] = start, shrinkedInterval[1] = end;

    // return
    return shrinkedInterval;
}

// Function to define error function 1
double firstErrorFunction(Polynomial polynomial, double *interval, double x) {
    // variable declaration
    double start = interval[0];
    double end = interval[1];
    double min; // the minimum value of |f'(x)| in [start, end]
    double startValue, endValue; 

    // get the first order derivative function
    Polynomial df = getDerivativeFunction(polynomial);

    // find min
    startValue = fabs(getValue(df, start));
    endValue = fabs(getValue(df, end));
    min = (startValue < endValue) ? startValue : endValue; 

    // compute and return the value of the error
    return fabs(getValue(polynomial, x)) / min;
}

// Function to define error function 2
double secondErrorFunction(Polynomial polynomial, double *interval, double cur, double prev) {
    // variable declaration
    double start = interval[0];
    double end = interval[1];
    double min; // the minimum value of |f'(x)| in [start, end]
    double max; // the maximum value of |f''(x)| in [start, end]

    // get the first order derivative function
    Polynomial df = getDerivativeFunction(polynomial);
    Polynomial ddf = getDerivativeFunction(df);

    // find min
    double firstOrderStartValue = fabs(getValue(df, start));
    double firstOrderEndValue = fabs(getValue(df, end));
    min = (firstOrderStartValue < firstOrderEndValue) ? firstOrderStartValue : firstOrderEndValue;

    // find max
    double secondOrderStartValue = fabs(getValue(ddf, start));
    double secondOrderEndValue = fabs(getValue(ddf, end));
    max = (secondOrderStartValue > secondOrderEndValue) ? secondOrderStartValue : secondOrderEndValue;

    // compute and return the value of the error
    return (max * pow(fabs(cur - prev), 2.0)) / (2 * min);
}

// Function to find a root with constant number of iteration
double findRootWithConstantIteration(Polynomial polynomial, double *interval) {
    // variable declaration
    double start = interval[0];
    double end = interval[1];
    int i; // loop variable
    int n; // number of iteration
    double x0 = (start + end) / 2;
    double x1; // next value for computation
    double firstError, secondError; // error for the first and second error function respectively

    // get the first order derivative of the polynomial
    Polynomial df = getDerivativeFunction(polynomial);

    // prompt user for number of iteration
    do {
        printf("Enter number of iteration: ");
        scanf("%d", &n);
        if (n <= 0) {
            printf("Number of iteration must be greater than 0!\n");
        }
    } while (n <= 0);

    // find root
    for (i = 0; i < n; i++) {
        double y = getValue(polynomial, x0);
        double dy = getValue(df, x0);

        // get next value 
        x1 = x0 - y/dy;

        // get the first error (with first error function) and second error (with second error function)
        firstError = firstErrorFunction(polynomial, interval, x1);
        secondError = secondErrorFunction(polynomial, interval, x1, x0);

        // store current value
        x0 = x1;
    }

    // display result
    printf("===== ROOT FOUND BY CONSTANT NUMBER OF ITERATION =====\n");
    printf("x = %lf\n\n", x0);
    printf("First error: %.20lf\n", firstError);
    printf("Second error: %.20lf\n", secondError);
    fprintf(output, "===== ROOT FOUND BY CONSTANT NUMBER OF ITERATION =====\n");
    fprintf(output, "Number of iteration: %d\n", n);
    fprintf(output, "x = %lf\n", x0);
    fprintf(output, "First error: %.20lf\n", firstError);
    fprintf(output, "Second error: %.20lf\n\n", secondError);

    return x0;
}

// Function to find a root with given error e and use the first error function
double findRootWithGivenErrorUseFirstErrorFunc(Polynomial polynomial, double *interval) {
    // variable declaration
    double e; // given error                       
    double start = interval[0];
    double end = interval[1];
    int i; // loop variable
    double x0 = (start + end) / 2;
    double error;

    // prompt user for error
    do
    {
        printf("Enter given error: ");
        scanf("%lf", &e);
        if (e <= 0)
        {
            printf("The given error must be a small number greater than 0!\n");
        }
    } while (e <= 0);

    // get the first order derivative of the polynomial
    Polynomial firstOrderDerivative = getDerivativeFunction(polynomial);

    // compute error using first error function
    error = firstErrorFunction(polynomial, interval, x0);

    // find root
    while (error > e) {
        double y = getValue(polynomial, x0);
        double dy = getValue(firstOrderDerivative, x0);

        // get next value
        x0 = x0 - y / dy;

        // computer error
        error = firstErrorFunction(polynomial, interval, x0);
    }

    // display result
    printf("===== ROOT FOUND BY FIRST ERROR FUNCTION =====\n");
    printf("x = %lf\n\n", x0);
    fprintf(output, "===== ROOT FOUND BY FIRST ERROR FUNCTION =====\n");
    fprintf(output, "Error: %lf\n", e);
    fprintf(output, "x = %.15lf\n\n", x0);

    // return 
    return x0;
}

// Function to find a root with given error e and use the second error function
double findRootWithGivenErrorUseSecondErrorFunc(Polynomial polynomial, double *interval) {
    // variable declaration
    double e; // given error
    double start = interval[0];
    double end = interval[1];
    int i; // loop variable
    double x0 = (start + end) / 2;
    double x1; // next value 
    double error; // error computed by second error function

    // prompt user for error
    do
    {
        printf("Enter given error: ");
        scanf("%lf", &e);
        if (e <= 0)
        {
            printf("The given error must be a small number greater than 0!\n");
        }
    } while (e <= 0);

    // get the first order derivative of the polynomial
    Polynomial firstOrderDerivative = getDerivativeFunction(polynomial);

    // find root
    i = 0;
    while (1)
    {
        double y = getValue(polynomial, x0);
        double dy = getValue(firstOrderDerivative, x0);

        // get next value
        x1 = x0 - y / dy;

        // computer error
        error = secondErrorFunction(polynomial, interval, x1, x0);

        // store current value
        x0 = x1;

        if (error < e) {
            break;
        }

        // counter
        i++;
    }

    // display result
    printf("===== ROOT FOUND BY SECOND ERROR FUNCTION =====\n");
    printf("x = %lf\n\n", x0);
    fprintf(output, "===== ROOT FOUND BY SECOND ERROR FUNCTION =====\n");
    fprintf(output, "Error: %lf\n", e);
    fprintf(output, "x = %.15lf\n\n", x0);

    // return
    return x0;
}

// Function to find a root with given error e satisfying |Xn - Xn-1| <= e
double findRootWithGivenError(Polynomial polynomial, double *interval) {
    // variable declaration
    double e; // given error
    double start = interval[0];
    double end = interval[1];
    int i; // loop variable
    double x0 = (start + end) / 2;
    double x1;    // next value

    // prompt user for error
    do
    {
        printf("Enter given error: ");
        scanf("%lf", &e);
        if (e <= 0)
        {
            printf("The given error must be a small number greater than 0!\n");
        }
    } while (e <= 0);

    // get the first order derivative of the polynomial
    Polynomial firstOrderDerivative = getDerivativeFunction(polynomial);

    // find root
    while (1)
    {
        double y = getValue(polynomial, x0);
        double dy = getValue(firstOrderDerivative, x0);

        // get next value
        x1 = x0 - y / dy;

        if (fabs(x1-x0) <= e) {
            break;
        }

        // store current value
        x0 = x1;
    }

    // display result
    printf("===== ROOT FOUND BY GIVEN ERROR =====\n");
    printf("x = %lf\n\n", x0);
    fprintf(output, "===== ROOT FOUND BY GIVEN ERROR =====\n");
    fprintf(output, "Error: %lf\n", e);
    fprintf(output, "x = %.15lf\n\n", x0);

    // return
    return x0;
}

// Function to check if an interval contains a real root
int hasRealRoot(Polynomial polynomial, double start, double end)
{
    return getValue(polynomial, start) * getValue(polynomial, end) < 0;
}

// Function to find upper bound of range of real roots
double findUpperBound(Polynomial polynomial) {
    // variable's declaration
    int i; // loop variable
    int k = -1; // index of the first negative coefficient
    double B; // maximum value of the absolute value among negative coefficient
    int A;    // maximum of absolute value of coefficients except the first coefficient

    if (polynomial.coefficient[0] > 0) {
        // find the first negative coefficient
        for (i = 0; i <= polynomial.degree; i++) {
            if (polynomial.coefficient[i] < 0) {
                k = i;
                break;
            }
        }
        // if there is no negative coefficient => upperbound = 0
        if (k == -1) {
            return 0;
        }
        // otherwise find B
        else {
            B = fabs(polynomial.coefficient[k]);
            for (i = k; i <= polynomial.degree; i++) {
                if (polynomial.coefficient[i] < 0 &&  B < fabs(polynomial.coefficient[i])) {
                    B = fabs(polynomial.coefficient[i]);
                }
            }
            return 1 + pow(fabs(B / polynomial.coefficient[0]), 1.0 / k);
        }
    }
    else {
        // find maximum of absolute value of coefficients except the first coefficient
        A = fabs(polynomial.coefficient[1]);
        for (i = 2; i < polynomial.degree; i++) {
            if (A < fabs(polynomial.coefficient[i])) {
                A = fabs(polynomial.coefficient[i]);
            }
        }

        return 1 + A / fabs(polynomial.coefficient[0]);
    }
}

// Function to find lower bound of range of real roots
double findLowerBound(Polynomial polynomial) {
    // variable's declaration
    int i;
    Polynomial polynomialTemp;
    polynomialTemp.coefficient = (double*)malloc(sizeof(double) * (polynomial.degree + 1));
    polynomialTemp.degree = polynomial.degree;

    // Find P(-x)
    if (polynomial.degree % 2 == 0) {
        for (i = 0; i < polynomial.degree; i++) {
            polynomialTemp.coefficient[i] = (i % 2 == 0) ? polynomial.coefficient[i] : -polynomial.coefficient[i];
        }
    }
    else {
        for (i = 0; i < polynomial.degree; i++){
            polynomialTemp.coefficient[i] = (i % 2 != 0) ? polynomial.coefficient[i] : -polynomial.coefficient[i];
        }
    }
    polynomialTemp.coefficient[polynomialTemp.degree] = polynomial.coefficient[polynomial.degree];

    // return lower bound
    return (-1) * findUpperBound(polynomialTemp);
}

// Gradient decent function
double gradientDescent(Polynomial df, double guess, double upperBound)
{
    // variable's declaration
    int sign;
    double leaningRate = 0.001;
    double eps = 1e-10;

    // find extrema
    if (getValue(df, guess) == 0) {
        return guess;
    }
    else if (getValue(df, guess) < 0) {
        sign = -1;
    }
    else {
        sign = 1;
    }

    while (fabs(getValue(df, guess)) > eps) {
        guess = guess + sign * leaningRate * getValue(df, guess);
        if (guess > upperBound) {
            return upperBound;
        }
    }
    return guess;
}

// find all interval containing roots
int findAllIntervalContainingRoots(Polynomial polynomial)
{
    // variables
    rootIntervals = (double*)malloc(sizeof(double) * 1);
    double upperBound = findUpperBound(polynomial);
    double lowerBound = findLowerBound(polynomial);
    rootIntervals[0] = lowerBound;
    int size = 1;
    double step = 1e-3;
    double i = lowerBound;
    int sign;

    // find derivative function
    Polynomial df = getDerivativeFunction(polynomial);

    // find extrema
    while (i < upperBound) {
        i = gradientDescent(df, i, upperBound);
        if (hasRealRoot(polynomial, rootIntervals[size-1], i)) {
            rootIntervals = (double*)realloc(rootIntervals, sizeof(double) * (size + 1));
            rootIntervals[size++] = i;
        }
        sign = (getValue(df, i) > 0) ? 1 : -1;
        do {
            i += step;
            if (i >= upperBound) {
                break;
            }
        } while ( (getValue(df, i) * sign > 0) && (fabs(getValue(df, i) > getValue(df, i + step))));
    }

    return size;
}

// Function to display interval
void displayInterval(int size)
{
    // check if polynomial has a root
    fprintf(output, "===== ALL INTERVALS CONTAINING ROOTS =====\n");
    if (size == 1) {
        printf("The function has no root\n\n");
        fprintf(output, "The function has no root\n\n");
    }
    else {
        printf("===== ALL INTERVALS CONTAINING ROOTS =====\n");
        for (int i = 0; i < size - 1; i++) {
            printf("[%lf, %lf]\n", rootIntervals[i], rootIntervals[i+1]);
            fprintf(output, "[%lf, %lf]\n", rootIntervals[i], rootIntervals[i + 1]);
        }
        fprintf(output, "\n");
    }
}

// Function to display shrinked interval
void displayShrinkedInterval(double* interval, int size)
{
    fprintf(output, "===== SHRINKED INTERVAL =====\n");
    for (int i = 0; i < size - 1; i++)
    {
        printf("===== SHRINKED INTERVAL =====\n");
        printf("[%lf, %lf]\n", interval[i], interval[i + 1]);
        fprintf(output, "[%lf, %lf]\n", interval[i], interval[i + 1]);
    }
    fprintf(output, "\n");
}