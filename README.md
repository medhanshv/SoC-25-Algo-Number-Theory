# SoC-25-Algo-Number-Theory
## Completed implementation of following functions:
#### def pair_gcd(a: int, b: int) -> int:

#### def pair_egcd(a: int, b: int) -> tuple[int, int, int]:

#### def gcd(*args: int) -> int:

#### def pair_lcm(a: int, b: int) -> int:

#### def lcm(*args: int) -> int:

#### def are_relatively_prime(a: int, b: int) -> bool:

#### def mod_inv(a: int, n: int) -> int:

#### def crt(a: list[int], n: list[int]) -> int:

#### def is_quadratic_residue_prime(a: int, p: int) -> int:

#### def is_quadratic_residue_prime_power(a: int, p: int, e: int) -> int:
#### floor_sqrt(x: int) -> int : 
#### is_perfect_power(x: int) -> bool :
#### is_prime(n: int) -> bool :
#### gen_prime(m : int) -> int : 
#### gen_k_bit_prime(k: int) -> int : 
#### factor(n: int) -> list[tuple[int, int]] :
#### euler_phi(n: int) -> int : 
#### Implemented a class QuotientPolynomialRing. 
#### It contains an instance variable called pi_generator which would be the the "quotienting polynomial", and an instance variable called element to represent the element of the ring.
#### __init__(self, poly: list[int], pi_gen: list[int]) -> None : 
#### A static method Add(poly1: QuotientPolynomialRing, poly2: QuotientPolynomialRing) -> QuotientPolynomialRing
#### A static method Sub(poly1: QuotientPolynomialRing, poly2: QuotientPolynomialRing) -> QuotientPolynomialRing  
#### A static method Mul(poly1: QuotientPolynomialRing, poly2: QuotientPolynomialRing) -> QuotientPolynomialRing  
#### A static method GCD(poly1: QuotientPolynomialRing, poly2: QuotientPolynomialRing) -> QuotientPolynomialRing 
#### A static method Inv(poly: QuotientPolynomialRing) -> QuotientPolynomialRing 
#### aks_test(n: int) -> bool : 
#### get_generator(p : int) -> int : 
#### discrete_log(x: int, g: int, p: int) -> int :
#### legendre_symbol(a: int, p: int) -> int: 
#### jacobi_symbol(a: int, n: int) -> int: 
#### modular_sqrt_prime(x: int, p: int) -> int : 
#### modular_sqrt_prime_power(x: int, p: int, e: int) -> int:
#### modular_sqrt(x: int, n: int) -> int: 
#### is_smooth(m: int, y: int) -> bool: 
#### probabilistic_dlog(x: int, g: int, p: int) -> int: 
#### probabilistic_factor(n: int) -> list[tuple[int, int]]:

completed relevant chapters from Victor Shoup
Created a video explaining basics of what goes into implementing these functions
