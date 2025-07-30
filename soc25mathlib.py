import random
from typing import Tuple,List

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
    a%=p
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
        
def gen_k_bit_prime(k: int) -> int:
    while True:
        n = random.getrandbits(k)
        msb=1<<(k-1)
        n=n|msb
        n=n|1
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
    def Add(poly1: "QuotientPolynomialRing", poly2: "QuotientPolynomialRing") -> "QuotientPolynomialRing":
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
    def Sub(poly1: "QuotientPolynomialRing", poly2: "QuotientPolynomialRing") -> "QuotientPolynomialRing":
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
    def Mul(poly1: "QuotientPolynomialRing", poly2: "QuotientPolynomialRing") -> "QuotientPolynomialRing":
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
    def GCD(poly1: "QuotientPolynomialRing", poly2: "QuotientPolynomialRing") -> "QuotientPolynomialRing":
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
    def poly_xgcd(f: list[int], g: list[int]) -> Tuple[list[int], list[int], list[int]]:
        def poly_rem(A: list[int], B: list[int]) -> list[int]:
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
            
 #hi nilabha, i have not corrected aks test    
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

def get_generator(p: int) -> int:
    if p == 2:
        raise ValueError("No generator exists for p = 2")
    phi= p-1
    for g in range(2, p):
        if not any(pow(g,phi//q,p) == 1 for q,a in factor(phi)):
            return g
    return -1 

def order(g: int ,p: int)-> int:
    r=p-1
    for q,a in factor(p-1):
        while r%q==0 and pow(g,r//q,p)==1:
            r//=q
    return r

def discrete_log(x: int, g: int, p: int) -> int:
    x%=p
    r=order(g,p)
    if pow(x,r,p)!=1:
        raise Exception("x^order is not congruent to 1, discrete log DNE")
    
    m=int(r**0.5)+1
    smol={}
    cur=1   
    for j in range(m):
        if cur not in smol:
            smol[cur] = j
        cur=(cur*g)%p

    ginv = pow(mod_inv(g,p),m,p)

# look for x*(g^{–m})^i in the smol table
    gamma=x
    for i in range(m):
        if gamma in smol:
            return i*m+smol[gamma]
        gamma=(gamma*ginv)%p

    raise Exception("Discrete log not found")

def legendre_symbol(a: int, p: int) -> int:
    a%=p
    if a==0:
        return 0
    return is_quadratic_residue_prime(a,p)

def jacobi_symbol(a: int, n: int) -> int:
    if n%2==0:
        raise Exception("n must be a positive odd integer")
    a%=n
    if a==0:
        return 0
    if n==1:
        return 1
    res=1
    for p,e in factor(n):
        leg=legendre_symbol(a, p)
        if leg==0:
            return 0
        # multiply in leg^e; since leg is +1 or -1, leg^e=leg iff e odd
        if e%2 == 1:
            res*=leg
    return res

def find_nonresidue(p: int)->int:
    z=2
    if(is_quadratic_residue_prime(z,p)!=1):
        return z
    while(True):    
        z=random.randint(2,p)
        if(is_quadratic_residue_prime(z,p)!=1):
            return z

def modular_sqrt_prime(x: int, p: int) -> int:
    if x%p==0:
        return 0
    if p==2:
        return x%p
    if(is_quadratic_residue_prime(x,p)!=1):
        raise Exception("square root does not exist")
    else:
        if(p%4==3):
            ans=pow(x,(p+1)//4,p)
            return min(ans,p-ans)
        # we need to apply tonelli shanks now(taken from wiki)
        Q=p-1
        s=0
        while(Q%2==0):
            s+=1
            Q//=2
        z=find_nonresidue(p)
        M=s
        c=pow(z,Q,p)
        t=pow(x,Q,p)
        R=pow(x,(Q+1)//2,p)
        while(t!=1):
            i=0
            temp=t
            while(temp!=1):
                i+=1
                temp=(temp*temp)%p
            pow2=2**(M-i-1)
            b=pow(c,pow2,p)
            M=i
            c=(b*b)%p
            t=(t*b*b)%p
            R=(R*b)%p
        ans=min(R,p-R) 
        return ans

def modular_sqrt_prime_power(x: int, p: int, e: int) -> int:
    if x%p**e==0:
        return 0
    if e==1:
        return modular_sqrt_prime(x, p)
    
    if is_quadratic_residue_prime_power(x,p,e) != 1:
            raise Exception("square root does not exist")
    
    if p==2:
        if e==2:
            return 1
        pe=1<<e
        x=x%pe
        root=-1
        for y in range(1,pe,2):
            if (y*y)%pe==x:
                root=y
                break
        return root

    r=modular_sqrt_prime(x%p,p)
    mod=p
    for i in range(1, e):
        nextmod=mod*p
        k=(r*r-x)//mod
        # solve 2r*t=-k (mod p)
        inv=mod_inv(2*r,p)
        t=(-k*inv)%p
        # lift
        r=r+t*mod
        mod=nextmod
    
    tmp = r%(p**e)
    ans=min(tmp,p**e-tmp)
    ans = int(ans)
    return ans

def modular_sqrt(x: int, n: int) -> int:
    x %= n
    if n==1:
        return 0
    facs=factor(n)  
    listmod= []
    mods= []
    for p,e in facs:
        pe = p**e
        if x%pe == 0:
            listmod.append([0])
        elif is_quadratic_residue_prime_power(x, p, e) == -1:
            raise Exception("No square root exists")
        else:
            r = modular_sqrt_prime_power(x, p, e)
            listmod.append([r,(-r) % pe])
        mods.append(pe)

    allroots = []
    def recurse(i: int, current: List[int]) -> int:
        if i == len(listmod):
            val = crt(current,mods) % n
            allroots.append(val)
            return 0
        for r in listmod[i]:
            recurse(i+1,current+[r])
        return 0
    
    recurse(0,[])
    return min(allroots)

def is_smooth(m: int, y: int) -> bool:
    if m<1:
        raise Exception("m must be a positive integer")
    if y<2:
        raise Exception("y must be at least 2")
    for p,e in factor(m):
        if p>y:
            return False
    return True

def probabilistic_factor(n: int) -> list[tuple[int,int]]:
    if n<=1:
        return []
    #our prime getter
    def primes_uptoy(y: int) -> list[int]:
        return [i for i in range(2,y+1) if is_prime(i)]
    #we factor m over a given list of (found)primes
    def div_exp(m: int, primes:list[int]) -> list[int]:
        exps=[]
        for p in primes:
            e=0
            while m%p == 0:
                m//=p
                e+=1
            exps.append(e)
        return exps
    #find F(2) dependency among vectors
    def gf2_depend(vecs: list[list[int]]) -> list[int]:
        N=len(vecs)
        basis: list[tuple[list[int], list[int]]] = []
        for i,v in enumerate(vecs):
            mask=[0]*N
            mask[i]=1
            w = v[:]
            for bvec,bmask in basis:
                pivot = next((j for j, bit in enumerate(bvec) if bit),None)
                if pivot is not None and w[pivot]:
                    w = [w[j]^bvec[j] for j in range(len(w))]
                    mask = [mask[j]^bmask[j] for j in range(N)]
            if all(bit == 0 for bit in w):
                return mask
            basis.append((w,mask))
        #we really should never get here if len(vecs)>dimension
        return []

    def sef_factor(m: int, y: int) -> int:
        sprimes = primes_uptoy(y)
        for p in sprimes:
            if m%p == 0:
                return p
        k = len(sprimes)
        while True:
            try:
                alphas,exponents,delta = [],[],None
                # collect k+2 smooth relations according to shoup
                for i in range(k+2):
                    while True:
                        a=random.randint(1,m-1)
                        d=random.randint(1,m-1) if i == 0 else delta
                        if d is not None:
                            rep = (pow(a,2,m)*d)%m or m
                        eis=div_exp(rep,sprimes)
                        alphas.append(a)
                        exponents.append(eis)
                        if i==0:
                            delta=d
                        break
                vs=[[e%2 for e in eis] + [1] for eis in exponents]
                c=gf2_depend(vs)
                # combine exponents
                combined = [0]*(k+1)
                for i,ci in enumerate(c):
                    if ci:
                        for j in range(k):
                            combined[j]+=exponents[i][j]
                        combined[k]+=1
                # build alpha,beta
                alpha=1
                for a,ci in zip(alphas,c):
                    if ci:
                        alpha = (alpha*a)%m
                beta=1
                for j in range(k):
                    beta=(beta*pow(sprimes[j],combined[j]//2,m))%m
                t=combined[k]//2 
                if t is None or delta is None or m is None:
                    continue
                deltap=pow(delta,t,m)
                invdelta=mod_inv(deltap,m)
                beta=(beta*invdelta)%m

                invbeta=mod_inv(beta,m)
                gamma=(alpha*invbeta)%m
                if gamma not in (1,m-1):
                    return gcd(gamma-1,m)
                #else retry
            except Exception:
                continue

    facs: dict[int,int] = {}
    remaining=n
    # we first strip factors of 2
    cnt = 0
    while (remaining&1) == 0:
        remaining>>=1
        cnt+=1
    if cnt:
        facs[2]=cnt

    # if the number is a perfect power
    if remaining>1 and is_perfect_power(remaining):
        # find base,exp such that base**exp==remaining
        maxb = remaining.bit_length()
        for b in range(2,maxb+1):
            a=int(round(remaining**(1.0/b)))
            if a>1 and pow(a,b)==remaining:
                base,exp = a,b
                break
        #we factor the base
        for p,e in probabilistic_factor(base):
            facs[p]=facs.get(p,0) + e*exp
        return sorted(facs.items())

    #prime case
    if remaining>1 and is_prime(remaining):
        facs[remaining] = facs.get(remaining,0)+1
        return sorted(facs.items())

    #extract via SEF and recurse
    f=sef_factor(remaining,100)
    g=remaining//f
    for p,e in probabilistic_factor(f):
        facs[p]=facs.get(p,0)+e
    for p,e in probabilistic_factor(g):
        facs[p]=facs.get(p,0)+e

    return sorted(facs.items())
