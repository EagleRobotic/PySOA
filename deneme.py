import numpy as np
from urdfpy import URDF
from IPython import embed

class Robot:
    def __init__(self, params):
        self.n = params[0]
        self.jdof = params[1]
        self.jtype = params[2]
        self.pgraph = params[3]

    def Rot(self,w,th):
        '''
        w: axis of rotation
        th: angle, radian
        '''
        R = 0
        w_ = np.zeros(shape=(3,3))
        if len(w)==3:
            w_[0,0:3] = np.array([0, -w[2], w[1]])
            w_[1,0:3] = np.array([w[2], 0, -w[0]])
            w_[2,0:3] = np.array([-w[1], w[0], 0])
            R = np.eye(3) + w_*np.sin(th) + w_*w_*(1-np.cos(th)) 
            #Rodriguez formula for the rotation matrix
        return R

    def tipBodies(self):
        p = self.pgraph
        ntips = 1 #Number of tip bodies
        tips = []
        for i in range(len(p)-1):
            if(p[i+1]<p[i]):
                tips.append(p[i])
                ntips+=1
        tips.append(np.max(p))
        k = np.max(p) #Number of bodies

        return tips, k

    def robotModel(self):
        x = np.zeros((3,sum(self.jdof)))
        x[0,:] = np.ones(sum(self.jdof))

        y = np.zeros((3,sum(self.jdof)))
        y[1,:] = np.ones(sum(self.jdof))

        z = np.zeros((3,sum(self.jdof)))
        z[2,:] = np.ones(sum(self.jdof))

        tips, k = self.tipBodies()
        

        ################################
        PHI = np.eye(6*k)
        PHI_t = np.zeros(shape=(6*len(tips), 6*k))
        PHI_td = np.zeros(shape=(6*len(tips), 6*k))
        FT = np.zeros(shape=(6*len(tips)))

        ################################
        V=VD=A=B=F = np.zeros(6*k)
        M = np.eye(6*k) #Mass matrix
        H = np.zeros(shape=(6*k, len(self.jdof)))

        #################################
        TH= np.zeros(shape=(np.sum(self.jdof)))
        THD= np.zeros(shape=(np.sum(self.jdof)))
        THDD= np.zeros(shape=(np.sum(self.jdof)))
        TAU = np.zeros(shape=(np.sum(self.jdof)))

        link_len = np.tile([1, 0, 0], (k+1,1))
        link_len = 0.2 * link_len
        lci = 0.5*link_len
        hi = np.tile([0, 0, 1], (k,1))
        h = np.zeros(shape=np.shape(hi))

        for i in range(len(self.jdof)):
            h[i][0:3] = hi[i][0]*x[:,i] + hi[i][1]*y[:,i] + hi[i][2]*z[:,i]

        return h, TH, link_len, x, y, z, lci

    def updateCoor(self):
        uc = [] #update coordinates
        tips, k = self.tipBodies()
        n1 = k + len(tips)

        params = self.robotModel()
        h = params[0]
        TH = params[1]
        li = params[2]
        lci = params[6]
        
        x = params[3]
        y = params[4]
        z = params[5]

        count = 0

        for i in reversed(range(0,k)):
            if(self.jdof[i]>0): #Is it seperation node or not?
                for c in reversed(range(self.jdof[i]-count+1,self.n-count)):
                    uc.append(c)
                    uca = np.asarray(uc)
                    if(self.jtype[c]==1):
                        R = self.Rot(h[c,0:3],TH[c])
                        x[:,uca] = np.dot(R,x[:,uca])
                        y[:,uca] = np.dot(R,y[:,uca])
                        z[:,uca] = np.dot(R,z[:,uca])
                    count +=1
                    if(self.jtype[c]==0):
                        ind = n1
                        li[ind,:] = li[ind,:] + h[c,:] * TH[c]
                n1 -=1
            if(self.jdof[i]==0):
                print("Seperation part!!!")


        link = np.array(np.shape(li))
        link_c = np.array(np.shape(lci))


        ii = 0
        cnt = 0
        for i in range(n1):
            if i==0:
                link[i,:] = link[i,:]
                link_c[i,:] = lci[i,:]
            else:
                ii = ii+self.jdof[cnt] #At the seperation node there is no increment
                if not self.jdof[cnt] == 0:
                    ind = ii
                    link[i,0:3] = li[i,0] * x[:,ind] + li[i,1] * y[:,ind] + li[i,2] * z[:,ind] 
                    link_c[i,:] = lci[i,0] * x[:,ind] + lci[i,1] *y[:,ind] + lci[i,2] * z[:,ind]
                    cnt +=1
                else:
                    #Jumps the seperation node and link update respect to previous node
                    jt = np.argwhere(self.pgraph==np.max(self.pgraph))
                    print("tryy")
        

    def readURDF(self):
        robot = URDF.load('my_robot.urdf')
        print(robot.links)


    #def create_phi(self):


def main():

    pgraph = np.array([1,2,3,4,5]) # Graph vector
    jdof = np.array([1,1,1,1,1]) # Joints degree of freedom
    #TODO: Design translational or rotational value assign later.
    jtype = np.array([1,1,1,1,1]) # For rotational joint:1, translation:0
    n = np.sum(jdof) #Total degree of freedom the body

    params = []
    params.append(n)
    params.append(jdof)
    params.append(jtype)
    params.append(pgraph)

    Robot1 = Robot(params)
    Robot.tipBodies(Robot1)
    Robot.robotModel(Robot1)
    Robot.updateCoor(Robot1)
    
    #TODO: It is not necassary for now
    #Robot.readURDF(Robot1)


if __name__ == "__main__" :
    main()