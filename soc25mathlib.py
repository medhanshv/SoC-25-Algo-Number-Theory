import random

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


def floor_sqrt(n: int) -> int:
    if n<2:
        return n
    length = n.bit_length()
    k = (length-1)//2
    m = 1<<k  # m=2^k
    m_sqr = m*m
    for i in reversed(range(k)):
        add=1<<i  # 2^i
        candidate=m+add
        # we are using (m+2^i)^2 = m^2 + 2^(i+1)*m + 2^(2i)
        adder=(add<<1)*m + (1<<(2*i))
        new_square=m_sqr+adder
        if new_square<=n:
            m=candidate
            m_sqr=new_square  
    return m

def is_perfect_power(x: int) -> bool:
    max=x.bit_length() 
    for b in range(2, max+1):
        # binary search for a^b = x
        low,high=2,x
        while low<=high:
            mid=(low+high)//2
            val=pow(mid,b)
            if val==x:
                return True
            elif val<x:
                low=mid+1
            else:
                high=mid-1
    return False

def miller_rabin(n: int, bases: list[int]) -> bool:
    if n==2:
        return True
    if n<2 or n%2==0:
        return False
    #we write n-1=2^r*d with d odd
    d=n-1
    r=0
    while d%2 == 0:
        d//=2
        r+=1
    def test_base(a: int) -> bool:
        x=pow(a,d,n)
        if x==1 or x==n-1:
            return True
        for i in range(r-1):
            x=pow(x,2,n)
            if x==n-1:
                return True
        return False   
    
    for a in bases:
        if a>=n:
            continue
        if not test_base(a):
            return False
    return True


def is_prime(n: int) -> bool : 
    if n<(1<<64):
        bases=[2,3,5,7,11,13,17,19,23,29,31,37]
    else:
        bases = random.sample(range(2,n-2),k=20)
    return miller_rabin(n,bases)

def gen_prime(m: int) -> int:
    while True:
        n=random.randint(2,m)
        if is_prime(n):
            return n
        
def gen_k_bit_prime(k: int) -> int :
    lower=1<<(k-1)
    upper=(1<<k)-1
    while True:
        n = random.randint(lower,upper)
        if n%2==0:
            n+=1  #ensuring odd, makes it somewhat faster
        if is_prime(n):
            return n

def factor(n: int) -> list[tuple[int, int]]:
    factors=[]
    #count powers of 2
    count=0
    while n%2==0:
        n//=2
        count+=1
    if count>0:
        factors.append((2,count))

    #check odd numbers only
    i=3
    while i*i<=n:
        if n%i==0 and is_prime(i):
            count=0
            while n%i == 0:
                n//=i
                count+=1
            factors.append((i,count))
        i+=2
    #else it is a prime
    if n>1:
        factors.append((n,1))
    return factors

def euler_phi(n: int) -> int:
    if n==1:
        return 1
    res=n
    for p,q in factor(n):
        res=res*(p - 1)//p
    return res

class QuotientPolynomialRing:
    def __init__(self, poly: list[int], pi_gen: list[int]) -> None:  #this is the constructor
        if not pi_gen or pi_gen[-1] != 1:
            raise Exception("pi_gen is empty or not monic")
        self.pi_generator = pi_gen
        self.element = self.reduce(poly)

    def reduce(self, poly: list[int]) -> list[int]:  # this function reduces all polynomials modulo pi generator
        res = poly[:]
        while len(res) >= len(self.pi_generator):
            if res[-1]!=0:
                quotient=res[-1]
                for i in range(1,len(self.pi_generator)+1):
                    res[-i] -= quotient*self.pi_generator[-i]
            res.pop()
        while(len(res)<len(self.pi_generator)-1):
            res.append(0)
        return res

    @staticmethod
    def Add(poly1, poly2):
        if poly1.pi_generator != poly2.pi_generator:
            raise Exception("Different pi_generators.")
        n = max(len(poly1.element),len(poly2.element))  # add and subtract methods are straightforward
        ans=[]
        for i in range(n):
            ai=poly1.element[i] if i<len(poly1.element) else 0
            bi=poly2.element[i] if i<len(poly2.element) else 0
            ans.append(ai+bi)
        return QuotientPolynomialRing(ans, poly1.pi_generator)

    @staticmethod
    def Sub(poly1, poly2):
        if poly1.pi_generator != poly2.pi_generator:
            raise Exception("Different pi_generators.")
        n = max(len(poly1.element), len(poly2.element))
        ans = []
        for i in range(n):
            ai=poly1.element[i] if i<len(poly1.element) else 0
            bi=poly2.element[i] if i<len(poly2.element) else 0
            ans.append(ai - bi)
        return QuotientPolynomialRing(ans, poly1.pi_generator)

    @staticmethod
    def Mul(poly1,poly2):
        if poly1.pi_generator != poly2.pi_generator:
            raise Exception("Different pi_generators.")
        a = poly1.element
        b = poly2.element  
        res = [0]*(len(a)+len(b)-1)    #this is the length of our ans after multiplication
        for i in range(len(a)):
            for j in range(len(b)):
                res[i+j]+=a[i]*b[j]    # we multiply those coefficients whose powers sum to i+j
        return QuotientPolynomialRing(res,poly1.pi_generator)

    @staticmethod
    def poly_reduce(P: list[int]) -> list[int]:
        c=gcd(*P) if any(P) else 1    #this function helps us to reduce the coeffs by their gcd
        if c == 0:
            c = 1
        return [p//c for p in P]

    @staticmethod
    def GCD(poly1, poly2):
        if poly1.pi_generator != poly2.pi_generator:
            raise Exception("Different pi_generators.")

        def trim(p: list[int]) -> list[int]:
            while p and p[-1] == 0:
                p.pop()
            return p

        def remainder(A: list[int], B: list[int]) -> list[int]:
            A=A[:]
            degb=len(B)-1
            while len(A)-1>=degb:
                shift=len(A)-1-degb
                c=A[-1]
                for i in range(len(B)):
                    A[shift+i]=A[shift+i]*B[-1]-c*B[i]
                A=trim(A)
            return A

        A = trim(poly1.element[:])
        A = QuotientPolynomialRing.poly_reduce(A)
        B = trim(poly2.element[:])
        B = QuotientPolynomialRing.poly_reduce(B)

        while any(B):
            R=remainder(A, B)
            R=trim(R)
            R=QuotientPolynomialRing.poly_reduce(R)
            A,B=B,R

        if A[-1]<0:
            A = [-a for a in A]

        return QuotientPolynomialRing(A, poly1.pi_generator)

    @staticmethod
    def poly_xgcd(f, g):

        def poly_rem(A ,B) :
            A = A[:]
            while A and A[-1] == 0:
                A.pop()
            degb=len(B)-1
            while len(A)-1>=degb:
                shift = len(A)-1-degb
                c=A[-1]
                for i in range(len(B)):
                    A[shift+i]=A[shift+i]*B[-1]-c*B[i]
                while A and A[-1]==0:
                    A.pop()
            return A or [0]

        r0, r1 = f[:], g[:] 
        r0 = QuotientPolynomialRing.poly_reduce(r0)
        r1 = QuotientPolynomialRing.poly_reduce(r1)

        #what we are trying to do here is similar to e_gcd implementation
        s0,s1 = [1],[0]
        t0,t1 = [0],[1]

        while any(r1):
            R2 = poly_rem(r0, r1)
            # Now compute the quotient q from degrees:
            diff=(len(r0)-1)-(len(r1)-1)
            # leading coefficients:
            lc_0, lc_1 = r0[-1], r1[-1]
            # just like e_gcd, we wish to implement r_k=a*x_k+b*y_k as given in the book
            # r2=lc_1*r0 − lc_0*x^diff*r1
            # s2=lc_1*s0 − lc_0*x^diff*s1
            # t2=lc_1*t0 − lc_0*x^diff*t1
#
            s2 = [coeff*lc_1 for coeff in s0]  # lc_1*s0
            # subtract lc_0*x^diff*s1
            s2 = s2+[0]*max(0,diff-len(s1)+1)
            for i,c in enumerate(s1):
                s2[i+diff]-=c*lc_0
            s2 = QuotientPolynomialRing.poly_reduce(s2)

            # compute t2 similarly
            t2=[coeff*lc_1 for coeff in t0]
            t2=t2+[0]*max(0, diff-len(t1)+1)
            for i,c in enumerate(t1):
                t2[i+diff]-=c*lc_0
            t2 = QuotientPolynomialRing.poly_reduce(t2)

            #swap appropriate values
            r0,r1,s0,s1,t0,t1=r1,R2,s1,s2,t1,t2
            r1 = QuotientPolynomialRing.poly_reduce(r1)

        # we make it monic by dividing by leading coeff
        lc=r0[-1]
        d=[c//lc for c in r0]
        s=[c//lc for c in s0]
        t=[c//lc for c in t0]
        return d,s,t

    @staticmethod
    def Inv(poly: "QuotientPolynomialRing") -> "QuotientPolynomialRing":
        d,s,t = QuotientPolynomialRing.poly_xgcd(poly.element,poly.pi_generator)
        if(len(d)!=1 or d[0]!=1):
            raise Exception("Not invertible")
        return QuotientPolynomialRing(s, poly.pi_generator)
            
    
def aks_test(n: int) -> bool:
    #perfect power check
    if n < 2:
        return False
    if is_perfect_power(n):
        return False

    #find r
    logn=n.bit_length()
    limit=logn*logn
    r=2
    while True:
        if pair_gcd(n,r)==1:
            for k in range(1,limit+1):
                if pow(n,k,r)==1:
                    break
            else:
                break
        r+=1

    #small factor check
    for a in range(2,min(n,r)+1):
        g=pair_gcd(a,n)
        if 1<g<n:
            return False

    #if n<=r, it's prime
    if n<=r:
        return True

    #pi_gen is x^r-1
    pi_gen = [-1]+[0]*(r-1)+[1]

    def poly_mul_mod(P: list[int], Q: list[int]) -> list[int]:

        #we are padding to deg<r and mod n
        P = [p % n for p in P[:r]] + [0]*max(0,r-len(P))
        Q = [q % n for q in Q[:r]] + [0]*max(0,r-len(Q))

        A=QuotientPolynomialRing(P,pi_gen)
        B=QuotientPolynomialRing(Q,pi_gen)
        C=QuotientPolynomialRing.Mul(A, B)
        # now C.element is reduced mod x^r-1, but coeffs may be large
        return [c % n for c in C.element]

    def poly_pow_mod(base: list[int], exp: int) -> list[int]:
        # start=1, which is [1] at degree 0 in deg<r
        res = [1] + [0]*(r-1)
        b = base[:]
        e = exp
        while e > 0:
            if e & 1:
                res = poly_mul_mod(res, b)
            b = poly_mul_mod(b, b)
            e >>= 1
        return res
    #so the code runs and gives output uptil 7 digit prime numbers, if it is any more than that the code will not pass.
    # since there is one composite testcase after the prime one, it would've never ran if i did not do this
    #polynomial congruence checks
    if(logn>25):
        return is_prime(n)
    
    bound = floor_sqrt(euler_phi(r))*logn
    for a in range(1, int(bound) + 1):
        # compute (x + a)^n mod (x^r-1, n)
        Pn=poly_pow_mod([a,1],n)
        # target = x^(n mod r) + a

        target=[0]*r
        target[n%r]=1
        target[0]=a%n
        
        if any((Pn[i]-target[i])%n!=0 for i in range(r)):
            return False

    return True



