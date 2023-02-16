import random
import time


FSS_SEC_PARA = 128
# FSS_SEC_PARA = 64
FSS_SEC_PARA_BYTE_LEN = int(FSS_SEC_PARA/8)

FSS_INPUT_LEN = 32
FSS_RING_LEN = 32
FSS_RING_BYTE_LEN = int(FSS_RING_LEN/8)


def sampleBits(seed,expand_len) -> int:
    if seed is None:
        random.seed(time.time())
    else:
        random.seed(seed)
    return random.getrandbits(expand_len)

class GroupElement(object):
    def __init__(self, value, bitlen,repr_value=None):
        assert (bitlen >= 1), "Improper bit length or scale"
       
        self.bitlen = bitlen
        self.Modulo = 2 ** self.bitlen
       
        if repr_value is None:
            self.value = (int(value) + 2 ** self.bitlen) % (2 ** self.bitlen)
        else:
            self.value = repr_value

    @classmethod
    def fromBin(cls,binary,bitlen):
        value = int.from_bytes(binary,"big")
        return GroupElement(value,bitlen)
    
    def getNegVal(self):
        return 2 ** self.bitlen - self.value


    def __add__(self, other):
        assert (type(other) is GroupElement), "Non groupType"
        assert (other.bitlen == self.bitlen), "can only be applied in the same bit length"

        value = (self.value + other.value) & (self.Modulo - 1)
        return GroupElement(value=None, bitlen=self.bitlen, repr_value=value)

    def __sub__(self, other):
        assert (type(other) is GroupElement), "Non groupType"
        assert (other.bitlen == self.bitlen), "can only be applied in the same bit length"

        value = (self.value - other.value+ self.Modulo) & (self.Modulo - 1)
        return GroupElement(value=None, bitlen=self.bitlen, repr_value=value)

    def __gt__(self, other):
        assert (type(other) is GroupElement), "Non groupType"
        assert (other.bitlen == self.bitlen), "can only be applied in the same bit length"
        return self.value > other.value
    
    def __lt__(self, other):
        assert (type(other) is GroupElement), "Non groupType"
        assert (other.bitlen == self.bitlen), "can only be applied in the same bit length"
        return self.value < other.value

    def __eq__(self, other):
        assert (type(other) is GroupElement), "Non groupType"
        assert (other.bitlen == self.bitlen), "can only be applied in the same bit length"
        return self.value == other.value

    def __getitem__(self, item):
        assert (self.bitlen >= item >= 0), f"No index at {item}"
        return self.value >> (self.bitlen-1-item) & 1
    
    def selfPrint(self):
        print("val is: ",self.value)
    
    def ele2Str(self):
        tmp=""
        for i in range(self.getLen()):
            tmp += str(self[i])
        return tmp
    
    def getLen(self):
        return self.bitlen

    def getValue(self):
        return self.value
    
    def packData(self):
        byteLen = int((self.bitlen+7)/8)
        # print("byteLen is: ",byteLen)
        return bytearray( self.value.to_bytes(byteLen,'big') )

class CW_DCF(object):
    def __init__(self, s, v_cw, t_l, t_r):
        self.s = s
        self.v_cw = v_cw
        self.t_l = t_l
        self.t_r = t_r
    
    def packData(self):
        binary = bytearray( self.s.to_bytes(FSS_SEC_PARA_BYTE_LEN,'big') )
        binary.extend( bytearray( self.v_cw.packData()  ) ) 
        binary.extend( bytearray( self.t_l.to_bytes(1,'big')  ) ) 
        binary.extend( bytearray( self.t_r.to_bytes(1,'big')  ) ) 
        return binary
    
    @classmethod
    def fromBin(cls,binary,ring_len):
        bytes_amount_per_cw = int((ring_len+7)/8)

        s = int.from_bytes(binary[:FSS_SEC_PARA_BYTE_LEN],'big')

        new_start = FSS_SEC_PARA_BYTE_LEN
        v_cw = GroupElement.fromBin(binary[new_start: new_start+bytes_amount_per_cw ], ring_len)
        new_start += bytes_amount_per_cw
        t_l = int.from_bytes(binary[new_start:new_start+1],'big')
        new_start += 1 
        t_r = int.from_bytes(binary[new_start:new_start+1],'big')
        return CW_DCF(s,v_cw,t_l,t_r)

class DCFKey(object):
    def __init__(self):
        self.seed = 0
        self.CW_List = []
        self.CW_payload = 0
    
    def packData(self):
        binary = bytearray( self.seed.to_bytes( FSS_SEC_PARA_BYTE_LEN,'big') )

        binary.extend( bytearray( self.CW_payload.packData() ))

        # Allow for at most 2**16 size 
        cw_size = len(self.CW_List)
        binary.extend( cw_size.to_bytes( 2,'big') )

        for v in self.CW_List:
            binary.extend( bytearray( v.packData()  ) ) 
        
        return binary
    
    @classmethod
    def fromBin(cls,binary,ring_len):
        dcfKey = DCFKey()

        bytes_amount_per_cw = int((ring_len+7)/8)


        dcfKey.seed = int.from_bytes(binary[:FSS_SEC_PARA_BYTE_LEN],'big')

        end = FSS_SEC_PARA_BYTE_LEN+bytes_amount_per_cw
        dcfKey.CW_payload = GroupElement.fromBin(binary[FSS_SEC_PARA_BYTE_LEN: end], ring_len)

        cw_size = int.from_bytes( binary[end:end+2 ],'big')

        new_start = end+2

        each_cw_len = int( len( binary[new_start: ] )/cw_size )

        for i in range(cw_size):
            start = new_start +i*each_cw_len
            end = new_start+(i+1)*each_cw_len
            dcfKey.CW_List.append(CW_DCF.fromBin( binary[start:end ],ring_len    ) )
        return dcfKey

class DDCFKey(object):
    def __init__(self):
        self.dcfKey =  DCFKey()
        self.Beta_share = 0

class DCF:
    """
    A DCF instantiation input output ring length, a given alpha value, output payload beta
    This functions returns DCF Key for if input < x, payload = 1 currently
    :param x:
    :param inverse: Keep False for unsigned comparison, keep None for signed comparison.
    :param sec_para:
    :param DEBUG:
    :return:
    """
    def __init__(self,sec_para=128,ring_len=32):
        """
        :param sec_para
        :param ring_len: Operation ring length
        :return:
        """
        self.sec_para = sec_para
        # Output ring length
        self.ring_len = ring_len

    def prg(self,seed):
        random.seed(seed)
        return random.getrandbits(4 * self.sec_para + 2)

    def convertG(self,_lambda):
        random.seed(_lambda)
        val =  random.getrandbits(self.ring_len)
        return GroupElement(val,self.ring_len)

    def keyGen(self,seed,alpha,beta) -> [DCFKey, DCFKey]:
        """
        This function returns DCF Key, where evaluation output to be "beta" if input < alpha, otherwise output 0
        :param bit_len: Operation group length
        :param alpha: A given alpha value with length may be not the same with output group element
        :param beta:  Desired payload value if input < alpha
        :return:
        """
        #seed_0 belongs to left, seed_1 belongs to right
        seed_bits = sampleBits(seed, self.sec_para*2)
        seed_0 = seed_bits >> self.sec_para
        seed_1 = seed_bits & ( (1<<self.sec_para) - 1)
        
        k0 = DCFKey()
        k1 = DCFKey()
        k0.seed = seed_0
        k1.seed = seed_1

        t_bits=[0,1]
        seeds=[seed_0,seed_1]
        
        V_a = GroupElement(0,self.ring_len)
        zero = GroupElement(0,self.ring_len)
        
        for i in range(alpha.getLen()):
            prg_res_l = self.prg(seeds[0])
            prg_res_r = self.prg(seeds[1])
            xor_res = prg_res_l ^ prg_res_r
            s_keep=[0,0]
            t_keep=[0,0]
            v_keep=[0,0]

            s_cw,v_cw=0,GroupElement(0,self.ring_len)
            # Line 10-12
            if alpha[i] == 0:
                s_cw = xor_res >> (self.sec_para+1) & ( (1<<self.sec_para) -1)
                v0_lose = prg_res_l >> 1 & ((1<<self.sec_para)-1)
                v1_lose = prg_res_r >> 1 & ((1<<self.sec_para)-1)
                v_cw = self.convertG(v1_lose) - self.convertG(v0_lose) - V_a
                if t_bits[1] == 1:
                    v_cw = zero - v_cw
                
                s_keep[0] = prg_res_l >> (3*self.sec_para+2)
                s_keep[1] = prg_res_r >> (3*self.sec_para+2)
                v_keep[0] = prg_res_l >> (2*self.sec_para+2) & ((1<<self.sec_para)-1)
                v_keep[1] = prg_res_r >> (2*self.sec_para+2) & ((1<<self.sec_para)-1)
                t_keep[0] = prg_res_l >> (2*self.sec_para+1) & 1
                t_keep[1] = prg_res_r >> (2*self.sec_para+1) & 1
            else:
                s_cw = xor_res >> (3*self.sec_para+2)
                v0_lose = prg_res_l >> (2*self.sec_para+2) & ((1<<self.sec_para)-1)
                v1_lose = prg_res_r >> (2*self.sec_para+2) & ((1<<self.sec_para)-1)
                v_cw = self.convertG(v1_lose) - self.convertG(v0_lose) - V_a
                if t_bits[1] == 1:
                    v_cw = zero - v_cw
                    v_cw -= beta
                else:
                    v_cw += beta
                
                s_keep[0] = prg_res_l >> (self.sec_para+1) & ((1<<self.sec_para)-1)
                s_keep[1] = prg_res_r >> (self.sec_para+1) & ((1<<self.sec_para)-1)
                v_keep[0] = prg_res_l >> 1 & ((1<<self.sec_para)-1)
                v_keep[1] = prg_res_r >> 1 & ((1<<self.sec_para)-1)
                t_keep[0] = prg_res_l & 1
                t_keep[1] = prg_res_r & 1
            
            # Line 14
            V_a -= self.convertG(v_keep[1])
            V_a += self.convertG(v_keep[0])
            #Line 14:last term addition for V_a
            if t_bits[1] == 1:
                V_a -= v_cw
            else:
                V_a += v_cw
            
            # Line 15-16
            t_cw_l = (xor_res >> (2*self.sec_para+1) & 1) ^ (1 - alpha[i])
            t_cw_r = (xor_res & 1) ^ alpha[i]
            CW = CW_DCF( s_cw, v_cw, t_cw_l, t_cw_r)
            t_cw_keep = t_cw_l
            if alpha[i] == 1:
                t_cw_keep = t_cw_r

            # Line 17-18
            for j in range(2):
                if t_bits[j] == 1:
                    seeds[j] = s_keep[j] ^ s_cw
                    t_bits[j] = t_keep[j] ^ t_cw_keep
                else:
                    seeds[j] = s_keep[j]
                    t_bits[j] = t_keep[j]
            k0.CW_List.append(CW)
            k1.CW_List.append(CW)

        payload = GroupElement(0,self.ring_len)
        if t_bits[1] == 0:
            payload += self.convertG(seeds[1])
            payload -= self.convertG(seeds[0])
            payload -= V_a
        else:
            payload -= self.convertG(seeds[1])
            payload += self.convertG(seeds[0])
            payload += V_a

        k0.CW_payload = payload
        k1.CW_payload = payload
        return k0, k1
    
    def eval(self, _id, x, key):
        """
        This function evaluates DCF at key with public value x
        :param party:
        :param x: a group element
        :param key:
        :return:
        """
        seed = key.seed
        t_bit = _id

        out = GroupElement(0,self.ring_len)

        levels = len(key.CW_List)
        for i in range(levels):
            # Line 3-4
            cw = key.CW_List[i]
            s_cw = cw.s
            v_cw = cw.v_cw
            t_cw_l = cw.t_l
            t_cw_r = cw.t_r
            prg_res = self.prg(seed)

            # Line 5-6
            s_l = prg_res >> (3*self.sec_para+2) 
            v_l = prg_res >> (2*self.sec_para+2) & ((1<<self.sec_para)-1)
            t_l = prg_res >> (2*self.sec_para+1) & 1

            s_r = prg_res >> (self.sec_para+1) & ((1<<self.sec_para)-1)
            v_r = prg_res >> 1 & ((1<<self.sec_para)-1)
            t_r = prg_res & 1

            if t_bit == 1:
                s_l ^= s_cw
                t_l ^= t_cw_l

                s_r ^= s_cw
                t_r ^= t_cw_r

            # Line 7-10
            if x[i] == 0:
                tmp = self.convertG(v_l)
                if t_bit==1:
                    tmp += v_cw
                if _id == 1:
                    out -= tmp
                else:
                    out += tmp
                seed = s_l
                t_bit = t_l
            else:
                tmp = self.convertG(v_r)
                if t_bit==1:
                    tmp += v_cw
                if _id == 1:
                    out -= tmp
                else:
                    out += tmp
                seed = s_r
                t_bit = t_r

        tmp = self.convertG(seed)
        if t_bit==1:
            tmp += key.CW_payload
        if _id == 1:
            out -= tmp
        else:
            out += tmp

        return out

class DDCF:
    """
    A DDCF instantiation input output ring length, a given alpha value, output payload beta
    This functions returns DCF Key for if input < x, payload = 1 currently
    :param x:
    :param inverse: Keep False for unsigned comparison, keep None for signed comparison.
    :param sec_para:
    :param DEBUG:
    :return:
    """
    def __init__(self,sec_para=128,ring_len=32):
        """
        :param sec_para
        :param ring_len: Operation ring length
        :return:
        """
        self.sec_para = sec_para
        # Output ring length
        self.ring_len = ring_len

    def convertG(self,_lambda):
        random.seed(_lambda)
        val =  random.getrandbits(self.ring_len)
        return GroupElement(val,self.ring_len)

    def keyGen(self,seed,alpha,beta1,beta2) -> [DDCFKey, DDCFKey]:
        """
        This function returns DCF Key, where evaluation output to be "beta" if input < alpha, otherwise output 0
        :param bit_len: Operation group length
        :param alpha: A given alpha value with length may be not the same with output group element
        :param beta:  Desired payload value if input < alpha
        :return:
        """
        beta = beta1 - beta2
        dcf = DCF( sec_para=self.sec_para, ring_len=self.ring_len)
        dcf_k0,dcf_k1 = dcf.keyGen(seed,alpha,beta)

        Beta_share0 = self.convertG(seed)
        Beta_share1 = beta2 - Beta_share0
        
        k0 = DDCFKey()
        k1 = DDCFKey()
        k0.dcfKey = dcf_k0
        k1.dcfKey = dcf_k1

        k0.Beta_share = Beta_share0
        k1.Beta_share = Beta_share1

        return k0,k1
    
    def eval(self, _id, x, key):
        """
        This function evaluates DCF at key with public value x
        :param party:
        :param x:
        :param key:
        :return:
        """
        dcf = DCF( sec_para=self.sec_para, ring_len=self.ring_len)
        dcfKey = key.dcfKey
        out = dcf.eval(_id,x,dcfKey)
        out += key.Beta_share
        return out

class BinIC:
    """
    Integer Comparison: Given an group element x, if 0<x<N/2 return 1, otherwise return 0, |x|< N/2.
    """
    def __init__(self,sec_para=128):
        self.sec_para = sec_para
        self.ddcf = DDCF(ring_len = 1)
    
    def keyGen(self,seed,input_len):
        self.input_len = input_len

        delta_val = sampleBits(seed, self.input_len)
        y_val = (1<<self.input_len) - delta_val

        # Get bits from pos [0..n-2]
        a = y_val & ((1<<(self.input_len-1)) -1)

        # Get bit by index n-1 Highest bit
        y_bit = (y_val>>(self.input_len-1)) & 1

        # print("y_bit: ", y_bit)
        # print("a: ",a)
        # print("y_val: ",y_val)

        delta = GroupElement(delta_val , self.input_len)
        # delta.selfPrint()

        d0 = GroupElement(sampleBits(None, self.input_len) , self.input_len)
        d1 = delta-d0


        # Convert a to a group element
        a_ele = GroupElement( a, self.input_len-1)

        beta1= GroupElement(1-y_bit,1)
        beta2= GroupElement(y_bit,1)

        k0, k1 = self.ddcf.keyGen(seed,a_ele,beta1,beta2)
        # zz = GroupElement( 2000000000, self.ring_len-1)
        # tmp = ddcf.eval(0,zz,k0)
        # tmp.selfPrint()

        # tmp += ddcf.eval(1,zz,k1)
        # tmp.selfPrint()
        return k0,d0,k1,d1
    
    # Start the online evaluation phase
    """
    param:: zeta is an group element
    param:: key is an BinIC key
    """
    def eval(self,_id, zeta, key):

        z_bit = zeta[0]
        z_len = zeta.getLen()

        zeta_low_bits = zeta.getValue() & ( (1<<(z_len-1)) -1)
        zz = GroupElement( (1<<(z_len-1)) - zeta_low_bits - 1, z_len-1)

        tmp1 = self.ddcf.eval(_id,zz,key)
        out = GroupElement(0,1)
        if _id==1:
            if z_bit==1:
                out = tmp1
            else:
                out = GroupElement( 1, 1) - tmp1
        else:
            if z_bit==1:
                out = tmp1
            else:
                out = GroupElement( 0, 1) - tmp1
        return out

class ICKey(object):
    """
    cw_payload
    keys: dcfKey
    """
    def __init__(self):
        self.CW_payload = 0
        self.keys = []
    
    def packData(self):
        binary = bytearray( self.CW_payload.packData() )
        for v in self.keys:
            binary.extend ( v.packData() )
        # print("ICKey len is: ",len(binary))
        return binary
    
    @classmethod
    def fromBin(cls,binary):
        icKey = ICKey()
        icKey.CW_payload = GroupElement.fromBin(binary[:FSS_RING_BYTE_LEN], FSS_RING_BYTE_LEN*8)
        dcfLen = int( len( binary[FSS_RING_BYTE_LEN: ])/2 )
        for i in range(2):
            start = FSS_RING_BYTE_LEN+i*dcfLen
            end = FSS_RING_BYTE_LEN+(i+1)*dcfLen
            icKey.keys.append(DCFKey.fromBin( binary[start:end ] ) )
        return icKey

class NewICKey(object):
    """
    cw_payload0: a arithmetical secret sharing of beta
    cw_payload1: correction word
    dcfKey
    """
    def __init__(self):
        self.CW_0 = 0
        self.CW_1 = 0
        self.dcf_key = DCFKey()
    
    def packData(self):
        binary = bytearray( self.CW_0.packData() )
        binary.extend ( self.CW_1.packData() )
        binary.extend ( self.dcf_key.packData() )
        return bytes( binary )
    
    @classmethod
    def fromBin(cls,binary,ring_len):
        # print("len binary is: ",len(binary))

        bytes_amount_per_cw = int((ring_len+7)/8)

        binary = bytearray(binary) #+ binary
        icKey = NewICKey()
        icKey.CW_0 = GroupElement.fromBin(binary[:bytes_amount_per_cw], ring_len)
        icKey.CW_1 = GroupElement.fromBin(binary[bytes_amount_per_cw:bytes_amount_per_cw*2], ring_len)
        icKey.dcf_key = DCFKey.fromBin( binary[bytes_amount_per_cw*2: ], ring_len )
        return icKey

class IntCmpNew:
    """
    Integer Comparison: 
        1. Interval test, if x \in [0,N/2] return b (b is an output group element), otherwise return 0.
        2. If keygen is given b1 and b2, goes similarly with one parameter b case, modify the payload value accordingly.
    """
    def __init__(self,sec_para=128,ring_len=32):
        self.sec_para = sec_para
        self.ring_len = ring_len
    
    def keyGen(self,seed,inputLen,b1,b2=None):
        beta = GroupElement(0, self.ring_len)
        addtional_payload = GroupElement(0, self.ring_len)

        if b2 is not None:
            beta = b1 - b2
            # addtional_payload = b2
        else:
            beta = b1

        #Calculate the first correction word
        cw0_0 = GroupElement(sampleBits(None, self.ring_len) , self.ring_len)
        cw0_1 = beta - cw0_0

        r_in = GroupElement(sampleBits(None, inputLen) , inputLen)
        r_in0 = GroupElement(sampleBits(None,inputLen) , inputLen)
        r_in1 = r_in - r_in0
        
        gamma = GroupElement( -1, inputLen)
        gamma += r_in

        dcf = DCF(sec_para=self.sec_para, ring_len = self.ring_len)
        dcfk0,dcfk1 = dcf.keyGen(seed,gamma,beta)

        #Calculate the second correction word
        alpha_p = r_in
        alpha_q = r_in + GroupElement( 1<<(inputLen-1), inputLen)
        alpha_q_prime = alpha_q + GroupElement( 1, inputLen) 

        scale = 0
        scale += 1 if alpha_p > alpha_q else 0
        scale -= 1 if alpha_p.getValue() > 0 else 0
        scale += 1 if alpha_q_prime.getValue() > ( (1<<(inputLen-1)) +1 ) else 0
        scale += 1 if alpha_q.getValue() == ( (1<<inputLen)-1 ) else 0

        # print("scale is: ",scale)
        scale *= beta.getValue()
        # To achieve general form of b1,b2 output
        cw_payload = GroupElement(scale, self.ring_len) + addtional_payload

        k0 = NewICKey()
        k0.CW_0 = cw0_0
        k0.dcf_key = dcfk0
        k0.CW_1 = GroupElement(sampleBits(None, self.ring_len) , self.ring_len)

        k1 = NewICKey()
        k1.CW_0 = cw0_1
        k1.dcf_key = dcfk1
        k1.CW_1 = cw_payload - k0.CW_1

        return r_in0,r_in1,k0,k1
    
    # Start the online evaluation phase
    """
    param:: zeta is an masked integer value
    param:: key is an IntCmpNew key
    """
    def eval(self,_id, zeta, key):
        dcf = DCF( sec_para=self.sec_para, ring_len = self.ring_len )

        inputLen = zeta.getLen()

        # start_ = time.time()
        scale = 1 if zeta.getValue() > 0 else 0

        # print("eval scale is: ",scale)
        scale -= 1 if zeta.getValue() > ( (1<<(inputLen-1)) +1 ) else 0

        # print("eval scale is: ",scale)

        scale *= key.CW_0.getValue()
        out = GroupElement(scale, self.ring_len)

        x_p = zeta + GroupElement( (1<<inputLen)-1, inputLen)
        x_q_prime = zeta + GroupElement( (1<<(inputLen-1))-2, inputLen)

        out -= dcf.eval(_id,x_p,key.dcf_key)
        out += dcf.eval(_id,x_q_prime,key.dcf_key)
        out += key.CW_1
        
        # print("One IntCmpNew eval cost ",time.time() - start_)

        return out

class IntCmp:
    """
    Integer Comparison: 
        1. Interval test, if x \in [0,N/2] return b (b is an output group element), otherwise return 0.
        2. If keygen is given b1 and b2, goes similarly with one parameter b case, modify the payload value accordingly.
    """
    def __init__(self,sec_para=128,ring_len=32):
        self.sec_para = sec_para
        self.ring_len = ring_len
    
    def keyGen(self,seed,inputLen, b1,b2=None):
        beta = GroupElement(0, self.ring_len)
        addtional_payload = GroupElement(0, self.ring_len)

        if b2 is not None:
            beta = b1 - b2
            addtional_payload = b2
        else:
            beta = b1

        seed_bits = sampleBits(seed, self.sec_para*2)
        seed_0 = seed_bits >> self.sec_para
        seed_1 = seed_bits & ( (1<<self.sec_para) - 1)

        dcf = DCF(sec_para=self.sec_para, ring_len = self.ring_len)
        k0 = ICKey()
        k1 = ICKey()

        r_in = GroupElement(sampleBits(None, inputLen), inputLen)
        r_in0 = GroupElement(sampleBits(None, inputLen) , inputLen)
        r_in1 = r_in - r_in0
        
        beta1= GroupElement(0 , self.ring_len) - beta
        dcfk0,dcfk1 = dcf.keyGen(seed_0,r_in,beta1)
        k0.keys.append(dcfk0)
        k1.keys.append(dcfk1)

        beta2 =  beta
        r_inEnd = GroupElement(1<<(inputLen-1), inputLen) + r_in
        dcfk00,dcfk11 = dcf.keyGen(seed_1,r_inEnd,beta2)
        k0.keys.append(dcfk00)
        k1.keys.append(dcfk11)

        cw_payload = GroupElement(0 , self.ring_len)

        if r_in > r_inEnd:
            cw_payload = GroupElement(beta.getValue() , self.ring_len)
        
        # To achieve general form of b1,b2 output
        cw_payload += addtional_payload
        
        k0.CW_payload = GroupElement(sampleBits(None, self.ring_len) , self.ring_len)
        k1.CW_payload = cw_payload - k0.CW_payload
        return r_in0,r_in1,k0,k1
    
    # Start the online evaluation phase
    """
    param:: zeta is an masked integer value
    param:: key is an IC key
    """
    def eval(self,_id, zeta, key):
        dcf = DCF( sec_para=self.sec_para, ring_len = self.ring_len )
        out = GroupElement( 0, self.ring_len)
        for i in range(2):
            out += dcf.eval(_id,zeta,key.keys[i])
        out += key.CW_payload
        return out

class FlexKey(object):
    """
    FlexKey defining what's transmited in OT.
    icKey: integer comparison key
    _id: evaluation id
    sign: 0 denotes evaluatinon result out*(-1), 1 denotes nothing
    payload: further addtional correction word, a group element
    """
    def __init__(self,icKey, _id ,sign, payload):
        self._id = _id
        self.sign = sign
        self.CW_payload = payload
        self.key = icKey

    """
    Constructing FlexKey from bytearray
    binary: a bytearray object
    """
    @classmethod
    def fromBin(cls,binary,ring_len):

        bytes_amount_per_cw = int((ring_len+7)/8)


        _id = int.from_bytes( binary[:1],"big" )
        sign = int.from_bytes( binary[1:2],"big" )
        CW_payload = GroupElement.fromBin(binary[2:2+bytes_amount_per_cw], ring_len)
        key = NewICKey.fromBin( binary[2+bytes_amount_per_cw:],ring_len )
        return FlexKey(key, _id, sign, CW_payload)
    
    def packData(self):
        """
            Byte length for each field
            _id  : 1
            sign : 1
            CW_payload : 
        """
        binary = bytearray( self._id.to_bytes(1,'big') )
        binary.extend( bytearray( self.sign.to_bytes(1,'big') ) )
        binary.extend( bytearray( self.CW_payload.packData()  ) ) 
        binary.extend ( self.key.packData() )
        # print("FlexKey len is: ",len(binary))
        return binary

class FlexFSS:
    """
    FlexFSS for (LOGICAL AND) computation with a IntCmp comparison result.
    Two types of FSS keys are provided for computing (A and B)
    0. Normal (A and B) keys
    1. Trap keys
    """
    def __init__(self,seed=None,sec_para=128,ring_len=32,fss_amount = 20):
        self.seed = seed
        self.sec_para = sec_para
        self.ring_len = ring_len
        self.refreshSeed()

        self.typesList=[]
        for i in range(fss_amount):
            new_t = random.randint(0,1)
            # print("fss t is: ",new_t)

            # Type 0 <-> normal keys
            # Type 1 <-> trap keys
            new_beta = GroupElement(sampleBits(None, self.ring_len) , self.ring_len)
            self.typesList.append( (new_t,new_beta) )

    def refreshSeed(self):
        seed_bits = sampleBits(self.seed, self.sec_para*2)
        self.seed = seed_bits & ( (1<<self.sec_para) - 1)
        random.seed(self.seed)

    """
    Prepare FSS keys for first evaluation, which is a combination of random normal/trap fss keys (FSS key and random output)
    """
    def initKeys(self):
        ic = IntCmpNew(sec_para=self.sec_para, ring_len = 1)
        initialKeys = []
        for t in self.typesList:
            # Refresh self.seed
            self.refreshSeed()

            # Type 0 <-> normal keys, set initial target output to be 1
            # Type 1 <-> trap keys, whatever the input, target output always be 0
            _type = t[0]

            # print("fss t is: ",t[0])

            beta = GroupElement(1, 1)
            if _type == 1:
                beta = GroupElement(0, 1)
            r0,r1,k0,k1 = ic.keyGen(self.seed,FSS_INPUT_LEN,beta)
            initialKeys.append( (r0,r1,k0,k1) )
        return initialKeys

    """
    Prepare message pairs for central bank
    """
    def keyGen4OT(self):
        ic = IntCmpNew(sec_para=self.sec_para, ring_len = self.ring_len)

        r_Array = []
        # Message pairs hold by two servers in OT transfer
        player0=[]
        player1=[]

        for t in self.typesList:
            # print("fss t is: ",t[0])
            self.refreshSeed()
            # First, preapare the two flexible fss keys
            # Type 0 <-> normal keys, set desired target output to be (beta,0)
            # Type 1 <-> trap keys, set initial target output to be (0,0)
            _type = t[0]
            beta = t[1] if _type == 0 else GroupElement( 0, self.ring_len )
            r0,r1,k0,k1 = ic.keyGen(self.seed,FSS_INPUT_LEN,beta)
            r_Array.append( r0+r1 )

            # Second, prepare OT messages
            add = t[1] if _type == 1 else GroupElement(0, self.ring_len)
            self.refreshSeed()
            payload0 = GroupElement( sampleBits(self.seed, self.ring_len), self.ring_len )
            payload1 = GroupElement(0, self.ring_len) - payload0

            flexKey0_0 = FlexKey(k0,0,1,payload0).packData()
            flexKey0_1 = FlexKey(k1,1,0,payload0 + add).packData()
            flexKey1_0 = FlexKey(k0,0,0,payload1 + add).packData()
            flexKey1_1 = FlexKey(k1,1,1,payload1).packData()

            
            player0.append( (flexKey0_0,flexKey0_1) )
            player1.append( (flexKey1_0,flexKey1_1) )
        return r_Array,player0,player1
    
    def getFSStypes(self):
        return  self.typesList
        # return [ele[0] for ele in self.typesList]

class FlexEval:
    """
    A virtual class to help evalute flex keys
    """
    def __init__(self,sec_para=128, ring_len = 32):
        self.sec_para = sec_para
        self.ring_len = ring_len
    
    def eval(self, flex: FlexKey, x):
        icKey = flex.key
        _id = flex._id
        sign = flex.sign
        payload = flex.CW_payload

        out = GroupElement( 0, self.ring_len)
        ic = IntCmpNew(sec_para=self.sec_para, ring_len = self.ring_len)
        if sign == 0:
            out -= ic.eval(_id, x, icKey)
        else:
            out += ic.eval(_id, x, icKey)
        out+= payload
        return out