#discrete math simulation
import random
import visual
from math import sin, cos, pi
import matplotlib.pyplot as plt
from time import sleep

plot_type = True
visualize = False
epsilon = 3


def gen_Watts_Strogatz(num_verts, K, Beta):
    edges = list()
    
    for i in range(num_verts):
        for j in range(num_verts):
            if (abs(i-j)<=(K/2) or abs(i-j)>(num_verts-(K/2)-1)) and i!=j and (i,j) not in edges and (j,i) not in edges:
                if i<j:
                    edges.append((i,j))
                else:
                    edges.append((j,i))
    new_edges = list()
    for vert in range(num_verts):
        for i,j in edges:
            if vert == i and i<j:
                if random.random() < Beta:
                    while 1:
                        terminus = random.randint(0,num_verts-1)
                        if terminus!=vert and (vert,terminus) not in new_edges and (terminus,vert) not in new_edges:
                            if vert<terminus:
                                new_edges.append((vert,terminus))
                            else:
                                new_edges.append((terminus,vert))
                            break
                else:
                    new_edges.append((i,j))
    return new_edges
    
    
def gen_Barabasi_Albert(num_verts, m, initial_path_length):
    #the graph is initialized by starting with a path
    if m > initial_path_length:
        raise ValueError, "m cannot be greater than initial path length"
    edges = [(i,i+1) for i in range(initial_path_length-1)]
    degrees = [0]*num_verts
    for i in range(num_verts):
        degrees[i]=len([edge for edge in edges if i in edge])
    for vertex in range(initial_path_length, num_verts):
        hat = list()
        for index in range(num_verts):
            hat += [index] * degrees[index]
        for new_edge in range(m):
            while True:
                potential_edge = random.choice(hat)
                if (vertex, potential_edge) not in edges and (potential_edge, vertex) not in edges:
                    edges.append((vertex, potential_edge))
                    degrees[vertex]+=1
                    degrees[potential_edge]+=1
                    break
                
    return edges
    
    
def gen_infection(verts,edges,size):
    new_verts = verts[:]
    neighborhoods = list()
    for i in range(len(verts)):
        connected_verts = [edge for edge in edges if i in edge]
        neighborhoods.append((len(connected_verts),random.random(),i))
    random.shuffle(neighborhoods)
    neighborhoods.sort(reverse=True)
    for entry in range(size):
        new_verts[neighborhoods[entry][2]]=1
    return new_verts
    
    
def spread(verts_in,edges,mortality_rate):
    num=len(verts_in)
    timeout=num
    finish_time=None
    verts=verts_in[:]
    if visualize:
        window = graph_vis(verts,edges)
    for time in range(timeout):
        healthy_neighborhood = [neighbor for neighbor in verts if neighbor==0]

        new_verts=verts[:]
        for i in range(num):
            if verts[i] == 0:
                connected_verts = [edge[0] for edge in edges if i in edge and i!=edge[0]]+[edge[1] for edge in edges if i in edge and i!=edge[1]]
                live_neighborhood = [neighbor for neighbor in connected_verts if verts[neighbor]!=3]
                infected_neighborhood = [neighbor for neighbor in connected_verts if verts[neighbor]==1]
                if len(infected_neighborhood) >= (len(live_neighborhood)/2.0):
                    new_verts[i] = 1
            elif verts[i] == 1:
                if random.random()<mortality_rate:
                    new_verts[i] = 3
                else:
                    new_verts[i] = 2
        verts=new_verts
        if visualize:
            window.update_color(verts, edges)
            sleep(0.5)
        if len(healthy_neighborhood)==0:
            finish_time=time
            break
    if visualize:   
        window.frame.visible = False
    return finish_time,verts
    
            
class graph_vis(object):
    def __init__(self, verts, edges):
        self.frame = visual.frame()
        self.susceptible = (0,1,0)
        self.infected = (1,0,0)
        self.recovered = (1,1,0)
        self.purple = (1,0,1)
        self.colors=[(0,1,0),(1,0,0),(1,1,0),(1,0,1)]
        self.v_verts, self.v_edges = self.gen_graph(verts, edges)
        
        
    def gen_graph(self, verts, edges):
        v_verts = list()
        v_edges = list()
        num = len(verts)
        r = (num*5)/(2*pi)
        self.position = list()
        for i in range(num):
            new_r = random.gauss(r,r/100)
            self.position.append((new_r*sin(2*(pi*i)/num),new_r*cos(2*(pi*i)/num),0))
        for i in range(num):
            v_verts.append(visual.sphere(frame=self.frame, pos=self.position[i],color=self.colors[verts[i]]))
            #if verts[i] == 0:#susceptible=green
            #    v_verts.append(visual.sphere(pos=position[i],color=(0,1,0)))
            #elif verts[i] == 1:#infected=red
            #    v_verts.append(visual.sphere(pos=position[i],color=(1,0,0)))
            #elif verts[i] == 2:#recovered=yellow
            #    v_verts.append(visual.sphere(pos=position[i],color=(1,1,0)))
            #elif verts[i] == 3:#dead=purple
            #    v_verts.append(visual.sphere(pos=position[i],color=(1,0,1)))
        for i,j in edges:
            v=(self.position[j][0]-self.position[i][0],self.position[j][1]-self.position[i][1],self.position[j][2]-self.position[i][2])
            v_edges.append(visual.cylinder(frame=self.frame, pos=self.position[i],axis=v,radius=0.25))
        return v_verts, v_edges
            
    def update_color(self, verts, edges):
        if len(verts) != len(self.v_verts) or len(edges) != len(self.v_edges):
            print verts
            print len(self.v_verts)
            print edges
            print len(self.v_edges)
            raise ValueError, "number of verts or edges have changed"
        for i in range(len(verts)):
            self.v_verts[i].color=self.colors[verts[i]]
        
    
def min_dynamo_BA():
    mortality_points=41
    iterations=30
    previous=0
    results=list()
    for n in range(2,7):
        temp_results = [0]*mortality_points
        for iter in range(iterations):
            mins = list()
            num=random.randint((n*10)-5,(n*10)+5)
            edges = gen_Barabasi_Albert(num, 4, 6)
            #edges = gen_Watts_Strogatz(num, 10, 0.5)
            for pre_mrate in range(mortality_points):
                mrate=pre_mrate/float(mortality_points-1)
                sizes=list()
                success_flag = 0
                fail_flag = 0
                for i in range(num):
                    if i%2==1:
                        infected_num = previous+((i+1)/2)
                    if i%2==0:
                        infected_num = previous-(i/2)
                    verts = [0]*num
                    verts = gen_infection(verts,edges,i)
                    ftime,verts_out=spread(verts,edges,mrate)
                    if ftime!=None:
                        sizes.append(i)
                        success_flag += 1
                    else:
                        fail_flag += 1
                    if success_flag>epsilon and fail_flag>epsilon:
                        break
    
                previous = min(sizes)
                mins.append(min(sizes)/float(num*5))
            #results.append(mins)
            print mins
            temp_results = [temp_results[m]+mins[m] for m in range(mortality_points)]
        if plot_type:
            x = [location/float(mortality_points-1) for location in range(0,mortality_points)]
            plt.plot(x,temp_results)
    if plot_type:
        plt.show()
    return results
        
        
def min_dynamo_WS():
    mortality_points=41
    iterations=30
    previous=0
    results=list()
    for n in range(2,7):
        temp_results = [0]*mortality_points
        for iter in range(iterations):
            mins = list()
            num=random.randint((n*10)-5,(n*10)+5)
            #edges = gen_Barabasi_Albert(num, 4, 6)
            edges = gen_Watts_Strogatz(num, 10, 0.5)
            for pre_mrate in range(mortality_points):
                mrate=pre_mrate/float(mortality_points-1)
                sizes=list()
                success_flag = 0
                fail_flag = 0
                for i in range(num):
                    if i%2==1:
                        infected_num = previous+((i+1)/2)
                    if i%2==0:
                        infected_num = previous-(i/2)
                    verts = [0]*num
                    verts = gen_infection(verts,edges,i)
                    ftime,verts_out=spread(verts,edges,mrate)
                    if ftime!=None:
                        sizes.append(i)
                        success_flag += 1
                    else:
                        fail_flag += 1
                    if success_flag>epsilon and fail_flag>epsilon:
                        break
    
                previous = min(sizes)
                mins.append(min(sizes)/float(num*5))
            #results.append(mins)
            print mins
            temp_results = [temp_results[m]+mins[m] for m in range(mortality_points)]
        if plot_type:
            x = [location/float(mortality_points-1) for location in range(0,mortality_points)]
            plt.plot(x,temp_results)
    if plot_type:
        plt.show()
    return results
    
    
def mortality_sim_WS(percent,infection_type):
    mortality_points=41
    iterations=30
    decades = 5
    results=list()
    for n in range(2,2+decades):
        temp_results = [0]*mortality_points
        for iter in range(iterations):
            mins = list()
            num=random.randint((n*10)-5,(n*10)+5)
            #edges = gen_Barabasi_Albert(num, 4, 6)
            edges = gen_Watts_Strogatz(num, 10, 0.5)
            for pre_mrate in range(mortality_points):
                mrate=pre_mrate/float(mortality_points-1)
                verts = [0]*num
                if infection_type==1:
                    verts = gen_infection(verts,edges,i)
                elif infection_type==0:
                    indecies = range(num)
                    random.shuffle(indecies)
                    for i in range(int(percent*num)):
                        verts[indecies[i]]=1
                ftime,verts_out=spread(verts,edges,mrate)
                percent_dead = len([vertex for vertex in verts_out if vertex==3])/float(num)
                mins.append(percent_dead/float(iterations))
            #results.append(mins)
            temp_results = [temp_results[m]+mins[m] for m in range(mortality_points)]
        print temp_results
        if plot_type:
            x = [location/float(mortality_points-1) for location in range(0,mortality_points)]
            plt.plot(x,temp_results)
    if plot_type:
        plt.show()
    return results

def mortality_sim_BA(percent,infection_type):
    mortality_points=41
    iterations=30
    decades = 5
    results=list()
    for n in range(2,2+decades):
        temp_results = [0]*mortality_points
        for iter in range(iterations):
            mins = list()
            num=random.randint((n*10)-5,(n*10)+5)
            edges = gen_Barabasi_Albert(num, 4, 6)
            #edges = gen_Watts_Strogatz(num, 10, 0.5)
            for pre_mrate in range(mortality_points):
                mrate=pre_mrate/float(mortality_points-1)
                verts = [0]*num
                if infection_type==1:
                    verts = gen_infection(verts,edges,i)
                elif infection_type==0:
                    indecies = range(num)
                    random.shuffle(indecies)
                    for i in range(int(percent*num)):
                        verts[indecies[i]]=1
                ftime,verts_out=spread(verts,edges,mrate)
                percent_dead = len([vertex for vertex in verts_out if vertex==3])/float(num)
                mins.append(percent_dead/float(iterations))
            #results.append(mins)
            temp_results = [temp_results[m]+mins[m] for m in range(mortality_points)]
        print temp_results
        if plot_type:
            x = [location/float(mortality_points-1) for location in range(0,mortality_points)]
            plt.plot(x,temp_results)
    if plot_type:
        plt.show()
    return results
    
def main():
    mrate=0.5
    num=20
    verts = [0]*num
    edges = gen_Barabasi_Albert(num, 10, 10)
    print "finished generation"
    infected_verts = gen_infection(verts,edges,num/4)
    print "finished infection"
    ftime,verts_out = spread(infected_verts,edges,mrate)
    print ftime
    #vis_graph(verts_out,edges)
    
    
if __name__=="__main__":
    mortality_sim_BA(0.3,0)