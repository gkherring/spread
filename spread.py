#discrete math simulation
import random
import visual
from math import sin, cos, pi
import matplotlib.pyplot as plt
from matplotlib import pylab
from time import sleep
from datetime import datetime as dt
from numpy import std, mean, sqrt
             
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
            
class simulate(object):
    def __init__(self, min_size, max_size, percentage, sim_type="mortality", infection_type="random", graph_type="BA"):
        self.min = min_size
        self.max = max_size
        self.infection_percent = percentage
        self.infection_type = infection_type
        self.graph_type = graph_type
        self.plot = True
        self.visualize = False
        self.epsilon = 3
        if sim_type=="mortality":
            self.mortality(percentage)
        elif sim_type=="dynamo":
            self.min_dynamo()
        elif sim_type=="fixed_mort":
            self.dynamo_spread(percentage)
        
    def gen_Watts_Strogatz(self, num_verts, K, Beta):
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
    
    def gen_Barabasi_Albert(self,num_verts, m, initial_path_length):
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
        
    def gen_infection_degreeranked(self,verts,edges,size):
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
    
    def gen_infection_random(self,verts,edges,size):
        new_verts = verts[:]
        indecies = range(len(verts))
        random.shuffle(indecies)
        for i in range(size):
            new_verts[indecies[i]]=1
        return new_verts
    
    def gen_infection_eigenranked(self,verts,edges,size):
        new_verts = verts[:] #Needs to be writen
        return new_verts
    
    def gen_infection_betweennessranked(self,verts,edges,size):
        new_verts = verts[:]#Needs to be writen
        return new_verts
    
    def spread(self, verts_in, edges, lethality_rate):
        num=len(verts_in)
        timeout=num
        finish_time=None
        verts=verts_in[:]
        if self.visualize:
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
                    if random.random()<lethality_rate:
                        new_verts[i] = 3
                    else:
                        new_verts[i] = 2
            verts=new_verts
            if self.visualize:
                window.update_color(verts, edges)
                sleep(0.5)
            if len(healthy_neighborhood)==0:
                finish_time=time
                break
        if self.visualize:   
            window.frame.visible = False
        return finish_time,verts
    
    def min_dynamo(self):
        lethality_points=51
        iterations=150
        previous=0
        results=list()
        temp_results = [list() for i in range(lethality_points)]
        for iter in range(iterations):
            mins = list()
            num=random.randint(self.min,self.max)
            if self.graph_type == "BA":
                edges = self.gen_Barabasi_Albert(num, 4, 6)
            else:
                edges = self.gen_Watts_Strogatz(num, 10, 0.5)
            for pre_mrate in range(lethality_points):
                mrate=pre_mrate/float(lethality_points-1)
                sizes=list()
                success_flag = 0
                fail_flag = 0
                for i in range(num):
                    if i%2==1:
                        infected_num = previous+((i+1)/2)
                    if i%2==0:
                        infected_num = previous-(i/2)
                    verts = [0]*num
                    if self.infection_type=="degree":
                        verts = self.gen_infection_degreeranked(verts,edges,i)
                    elif self.infection_type=="random":
                        verts = self.gen_infection_random(verts,edges,i)
                    elif self.infection_type=="eigen":
                        verts = self.gen_infection_eigenranked(verts,edges,i)
                    elif self.infection_type=="betweenness":
                        verts = self.gen_infection_betweennessranked(verts,edges,i)
                    ftime,verts_out=self.spread(verts,edges,mrate)
                    if ftime!=None:
                        sizes.append(i)
                        success_flag += 1
                    else:
                        fail_flag += 1
                    if success_flag>self.epsilon and fail_flag>self.epsilon:
                        break
                previous = min(sizes)
                temp_results[pre_mrate].append(min(sizes))
            #results.append(mins)
            #print temp_results
            temp_results = [temp_results[m].append(mins[m]) for m in range(lethality_points)]
        if self.plot:
            x = [location/float(lethality_points-1) for location in range(0,lethality_points)]
            figure = plt.figure(1)
            ax1 = fig.add_axes((.1,.4,.8,.5))
            plt.clf()
            y = [mean(item) for item in temp_results]
            error = [2.58*std(item)/sqrt(len(item)) for item in temp_results]
            print y
            print error
            pylab.errorbar(x,y,yerr=error,fmt="ro")
            txt=", ".join(["Averaged over "+str(iterations)+" distinct graphs","99% Confidence Interval","Min Graph Size: "+str(self.min),"Max Graph Size: "+str(self.max)+"\n", "Infecting type: "+self.infection_type,"Graph type: "+self.graph_type,"Showing min-dynamo ceiling vs Lethality"])
            plt.title(txt, fontsize=10)
            plt.xlabel("Lethality")
            plt.ylabel("Percent of Graph")
            figure.savefig(", ".join([str(dt.now()),self.infection_type,self.graph_type])+"_min_dynamo_ceiling.png")
            #plt.show()
        return results
    
    def mortality(self,percent):
        lethality_points=51
        iterations=200
        results=list()
        percent_dead_results = [list() for i in range(lethality_points)]
        percent_infected_results = [list() for i in range(lethality_points)]
        for iter in range(iterations):
            dead = list()
            infected = list()
            num=random.randint(self.min,self.max)
            if self.graph_type == "BA":
                edges = self.gen_Barabasi_Albert(num, 4, 6)
            else:
                edges = self.gen_Watts_Strogatz(num, 10, 0.5)
            for pre_mrate in range(lethality_points):
                mrate=pre_mrate/float(lethality_points-1)
                verts = [0]*num
                if self.infection_type=="degree":
                    verts = self.gen_infection_degreeranked(verts,edges,int(percent*num))
                elif self.infection_type=="random":
                    verts = self.gen_infection_random(verts,edges,int(percent*num))
                elif self.infection_type=="eigen":
                    verts = self.gen_infection_eigenranked(verts,edges,int(percent*num))
                elif self.infection_type=="betweenness":
                    verts = self.gen_infection_betweennessranked(verts,edges,int(percent*num))
                ftime,verts_out=self.spread(verts,edges,mrate)
                percent_dead = len([vertex for vertex in verts_out if vertex==3])/float(num)
                percent_affect = len([vertex for vertex in verts_out if vertex!=0])/float(num)
                percent_dead_results[pre_mrate].append(percent_dead)
                percent_infected_results[pre_mrate].append(percent_affect)
            #percent_dead_results = [percent_dead_results[m].append(dead[m]) for m in range(lethality_points)]
            #percent_infected_results = [percent_infected_results[m].append(infected[m]) for m in range(lethality_points)]
        #print percent_dead_results
        #print percent_infected_results
        if self.plot:
            x = [location/float(lethality_points-1) for location in range(0,lethality_points)]
            figure = plt.figure(1)
            plt.clf()
            yone = [mean(item) for item in percent_dead_results]
            ytwo = [mean(item) for item in percent_infected_results]
            errorone = [2.58*std(item)/sqrt(len(item)) for item in percent_dead_results]
            errortwo = [2.58*std(item)/sqrt(len(item)) for item in percent_infected_results]
            one = pylab.errorbar(x,yone,yerr=errorone,fmt='ro')
            two = pylab.errorbar(x,ytwo,yerr=errortwo,fmt='bo')
            print yone, ytwo, errorone, errortwo
            txt = ", ".join(["Percentage of the graph initially infected: "+str(percent),"Averaged over "+str(iterations)+" distinct graphs\n","99% Confidence Interval","Min Graph Size: "+str(self.min),"Max Graph Size: "+str(self.max), "Infecting type: "+self.infection_type+"\n","Graph type: "+self.graph_type,"Mortality and Spread vs Lethality"])
            plt.title(txt, fontsize=10)
            plt.legend([one,two],["Dead","Infected"])
            plt.xlabel("Lethality")
            plt.ylabel("Percent of Graph")
            figure.savefig(", ".join([str(dt.now()),self.infection_type,self.graph_type])+"_mortality.png")
            #plt.show()
        return percent_dead_results, percent_infected_results

    def dynamo_spread(self,lethality_percent):
        infected_percentage_points=51
        iterations=200
        results=list()
        percent_dead_results = [[] for i in range(infected_percentage_points)]
        percent_infected_results = [[] for i in range(infected_percentage_points)]
        for iter in range(iterations):
            dead = list()
            infected = list()
            num=random.randint(self.min,self.max)
            if self.graph_type == "BA":
                edges = self.gen_Barabasi_Albert(num, 4, 6)
            else:
                edges = self.gen_Watts_Strogatz(num, 10, 0.5)
            for pre_infection_percent in range(infected_percentage_points):
                percent=pre_infection_percent/float(infected_percentage_points-1)
                verts = [0]*num
                if self.infection_type=="degree":
                    verts = self.gen_infection_degreeranked(verts,edges,int(percent*num))
                elif self.infection_type=="random":
                    verts = self.gen_infection_random(verts,edges,int(percent*num))
                elif self.infection_type=="eigen":
                    verts = self.gen_infection_eigenranked(verts,edges,int(percent*num))
                elif self.infection_type=="betweenness":
                    verts = self.gen_infection_betweennessranked(verts,edges,int(percent*num))
                ftime,verts_out=self.spread(verts,edges,lethality_percent)
                percent_dead = len([vertex for vertex in verts_out if vertex==3])/float(num)
                percent_affect = len([vertex for vertex in verts_out if vertex!=0])/float(num)
                percent_dead_results[pre_infection_percent].append(percent_dead)
                percent_infected_results[pre_infection_percent].append(percent_affect)
            #percent_dead_results = [percent_dead_results[m].append(dead[m]) for m in range(infected_percentage_points)]
            #percent_infected_results = [percent_infected_results[m].append(infected[m]) for m in range(infected_percentage_points)]
        #print percent_dead_results
        #print percent_infected_results
        if self.plot:
            x = [location/float(infected_percentage_points-1) for location in range(0,infected_percentage_points)]
            figure = plt.figure(1)
            ax1 = figure.add_axes((.1,.4,.8,.5))
            plt.clf()
            yone = [mean(item) for item in percent_dead_results]
            ytwo = [mean(item) for item in percent_infected_results]
            print percent_infected_results
            errorone = [2.58*std(item)/sqrt(len(item)) for item in percent_dead_results]
            errortwo = [2.58*std(item)/sqrt(len(item)) for item in percent_infected_results]
            print x
            print ytwo, errortwo
            one = pylab.errorbar(x,yone,yerr=errorone,fmt='ro')
            two = pylab.errorbar(x,ytwo,yerr=errortwo,fmt='bo')
            txt=", ".join(["Lethality of Infection: "+str(percent),"Averaged over "+str(iterations)+" distinct graphs\n","99% Confidence Interval","Min Graph Size: "+str(self.min),"Max Graph Size: "+str(self.max), "Infecting type: "+self.infection_type+"\n","Graph type: "+self.graph_type,"Mortality and Spread vs Infected Size"])
            plt.title(txt,fontsize=10)
            plt.legend([one,two],["Dead","Infected"])
            plt.xlabel("Percent Infected")
            plt.ylabel("Percent of Graph")
            figure.savefig(", ".join([str(dt.now()),self.infection_type,self.graph_type])+"_dynamo_spread.png")
            #plt.show()
        return percent_dead_results, percent_infected_results
    
def main():
    for i in range(11):
        simulate(20, 70, i/float(10),"fixed_mort", "degree", "BA")
    for i in range(11):
        simulate(20, 70, i/float(10),"mortality", "degree", "BA")
    for i in range(11):
        simulate(20, 70, i/float(10),"fixed_mort", "random", "BA")
    for i in range(11):
        simulate(20, 70, i/float(10),"mortality", "random", "BA")
    
    
if __name__=="__main__":
    main()
