# if __name__ == "main":
RING_LEN = 1
seed=135124

def dcf_test():
    from fss import DCF, GroupElement

    dcf = DCF(ring_len = RING_LEN)

    a = GroupElement(2**31+11,RING_LEN)
    b = GroupElement(3,RING_LEN)

    k0, k1 = dcf.keyGen(seed,a,b)
    x= GroupElement(2**31-13,RING_LEN)
    v0 = dcf.eval(0,x,k0)
    v1 = dcf.eval(1,x,k1)

    v0.selfPrint()
    v1.selfPrint()

    v0+=v1
    v0.selfPrint()


def ddcf_test():
    from fss import DDCF, GroupElement

    ddcf = DDCF(ring_len = RING_LEN)

    beta1 = GroupElement( 0, RING_LEN)
    beta2 = GroupElement( 3, RING_LEN)

    a = GroupElement(2**31+11,RING_LEN-1)
    # print(a.getLen())
    # print(ele2Str(a))

    k0, k1 = ddcf.keyGen(seed, a, beta1, beta2)

    # 2**31-1 = 2147483647
    # weird number = 2087559245
    # for x in range()

    x= GroupElement(2**31+10,RING_LEN-1)
    # print(ele2Str(x))
    v0 = ddcf.eval(0,x,k0)
    v1 = ddcf.eval(1,x,k1)

    v0.selfPrint()
    v1.selfPrint()

    v0+=v1
    v0.selfPrint()


def BinIC_test(val):
    from fss import BinIC, GroupElement

    ic = BinIC()

    input_Len = 32

    k0,d0,k1,d1 = ic.keyGen(seed,input_Len)
    zeta = GroupElement( val, input_Len)
    zeta += d0
    zeta += d1
    
    # z_val = zeta.getValue()

    # print("z_val is: ", z_val)
    v0 = ic.eval(0,zeta,k0)
    v1 = ic.eval(1,zeta,k1)

    v0.selfPrint()
    v1.selfPrint()

    v0+=v1
    v0.selfPrint()


def IC_test():
    from fss import IntCmp, GroupElement

    ic = IntCmp(ring_len = RING_LEN)

    inputLen = 32
    beta1 = GroupElement( 1, RING_LEN)
    beta2 = GroupElement( 6, RING_LEN)
    r0,r1,k0,k1 = ic.keyGen(seed,inputLen,beta1)
    # r0,r1,k0,k1 = ic.keyGen(seed,beta1,beta2)

    zeta = GroupElement( -10000, inputLen)
    # zeta = GroupElement( 10000, inputLen)
    zeta += r0
    zeta += r1
    
    v0 = ic.eval(0,zeta,k0)
    v1 = ic.eval(1,zeta,k1)

    # v0.selfPrint()
    # v1.selfPrint()

    v0+=v1
    v0.selfPrint()


def NewIC_test(val,betaInt):
    from fss import IntCmpNew, GroupElement

    inputLen = 32

    ic = IntCmpNew(ring_len = RING_LEN)
    beta1 = GroupElement( betaInt, RING_LEN)
    # beta2 = GroupElement( 51, RING_LEN)
    # r0,r1,k0,k1 = ic.keyGen(seed,beta1)
    r0,r1,k0,k1 = ic.keyGen(seed,inputLen,beta1)

    # zeta = GroupElement( -10000, RING_LEN)
    zeta = GroupElement( val, inputLen)
    zeta += r0
    zeta += r1
    
    v0 = ic.eval(0,zeta,k0)
    v1 = ic.eval(1,zeta,k1)

    # v0.selfPrint()
    # v1.selfPrint()

    v0+=v1
    return v0.getValue()

# dcf_test()
# ddcf_test()

# for i in range(100):
#     BinIC_test(-i)


# IC_test()

import random
for i in range(1):
    val = random.randint(0, 1<<32)

    realCmp = 0
    if val <= ( 1<<31):
        # print("positive",val)
        realCmp = 1
    else:
        # print("negative",(1<<32) -val )
        realCmp = 0
    
    beta = 1
    cmpRet = NewIC_test( val, beta )

    # if cmpRet!=0 and cmpRet!= 1 and cmpRet!= 3:
    # if cmpRet== 3*beta or cmpRet==0:
    #     print(cmpRet)

    # if cmpRet== 0:
        # print(cmpRet)

    # cmpRet %= 2
    if realCmp*beta != cmpRet:
        print(i,"cmp result not desired", "real: ",realCmp, "exp: ",cmpRet)