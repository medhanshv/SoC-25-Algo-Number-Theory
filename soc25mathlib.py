def pair_gcd(a: int, b: int) -> int:
    m,n=a,b    
    while n!=0:
        m=m%n   #essentially applying euclidean algorithm
        m,n=n,m
    return m

def pair_egcd(a: int, b: int) -> tuple[int, int, int]:
    m,n=a,b
    x,i=1,0
    y,j=0,1
    while n!=0:     #here we store i and j which satisfy the recurrence r_k=a*x_k+b*y_k at all steps.
        q=m//n      #when the algorithm terminates, r_k=gcd and x,y store the desired coefficients
        m=m%n
        m,x,y,n,i,j=n,i,j,m,x-i*q,y-j*q
    return x,y,m

def gcd(*args: int) -> int:
    res=args[0]
    for num in args[1:]:
        res=pair_gcd(res,num)
    return res

def pair_lcm(a: int, b: int) -> int:
    res=(a*b)//pair_gcd(a,b)
    return res

def lcm(*args: int) -> int:
    res=args[0]
    for num in args[1:]:
        res=pair_lcm(res,num)
    return res

def are_relatively_prime(a: int, b: int) -> bool:
    if(pair_gcd(a,b)==1):
        return True
    return False

def mod_inv(a: int, n: int) -> int:
    x,y,d=pair_egcd(a, n)      #here ax+by=gcd(a,b)=1 in this case, so x%n is our modular inverse
    if d!=1:
        raise Exception("a and n are not coprime")
    return x % n  

def crt(a: list[int], n: list[int]) -> int:
    prod=1
    x=0
    for num in n:
        prod*=num
    for i in range(len(n)):
        b=(prod//n[i])
        x+=b*mod_inv(b,n[i])*a[i]
    return x%prod

def is_quadratic_residue_prime(a: int, p: int) -> int:
    if(p==2): 
        return 1
    euler=pow(a,(p-1)//2,p)
    if(euler==(p-1)): 
        return -1
    return euler

def is_quadratic_residue_prime_power(a: int, p: int, e: int) -> int:
    if(pair_gcd(a,p)!=1):
        return 0
    
    if(p==2):
        if(e==1):
            return 1
        elif(e==2):
            return 1 if (a%4==1) else -1
        else:
            return 1 if (a%8==1) else -1
           
    return is_quadratic_residue_prime(a,p)

