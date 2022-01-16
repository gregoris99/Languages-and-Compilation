% --> TAKES APROXIMATELY 1.5 MINS TO FINISH <--
%load data
load('xy.mat')

% STEP 1: Generate random population of pop_size chromosomes:
%make matrix to hold populations
pop_size=150;
population = zeros(pop_size,101);

%populate population
current_pop=0;

%Make first chromosome of population by randomising the order of the elements of array cities
%population(1,:)=cities(randperm(numel(cities)));
population(1,1:100)=randperm(100,100);
current_pop=current_pop+1;

%populate the rest of the populatio
while ( current_pop < pop_size )
    
    %Randomise the order of the elements of array cities
    temp_chrom=randperm(100,100);
    
    current_pop=current_pop+1;
    population(current_pop,1:100)=temp_chrom;

end

% Step 2, Calculate fitnes of each chromosome in innitial population
for i=1:pop_size
    population(i,101)=calc_fit(population(i,1:100),xy);
end

%Start Loop
%decide generations
generations=4000;
cross_probability=0.2;
mutation_probability=0.1;

%Loop
for gen=1:generations

    %Loop population_size/2, to make new generation
    new_gen=zeros(pop_size,101);
    current_new_gen=0;
    for pop=1:(pop_size/2)
        %Step 3, Selection of 2 chromosomes from innitial population
        parents=select_turnmnt(pop_size, population(:,101));
        %parents=select_rulette(pop_size, population(:,101));
        %crossover probability
        r=rand(1);
        temp_children=zeros(2,100);
        if(r<cross_probability)
             %Step 4, Crossover
             temp_children(:,:)=Order_crossover(population(parents(1),:), population(parents(2),:) );
             %temp_children=PartMap_crossover(population(parents(1),:), population(parents(2),:));
        else
             temp_children(1,:)=population(parents(1),1:100);
             temp_children(2,:)=population(parents(2),1:100);
        end
        
        %mutation probability
        r1=rand(1);
        temp_childrenMute=zeros(2,100);
        if(r1<mutation_probability)
             %Step 5, Mutaion
             %temp_childrenMute=flipMutation(population(parents(1),:), population(parents(2),:));
             temp_childrenMute(:,:)=swapMutation(population(parents(1),1:100), population(parents(2),1:100) );   
        else
            temp_childrenMute=temp_children;
        end
        
        %add childer to new generation
        current_new_gen=current_new_gen+1;
        new_gen(current_new_gen,1:100)=temp_childrenMute(1,1:100);
        current_new_gen=current_new_gen+1;
        new_gen((current_new_gen),1:100)=temp_childrenMute(2,1:100);
    end
   
    %replace population with Elitism
    populationSorted = sortrows(population,101);
    %old_population_temp=population;
    population=[];
    population(1:20,:)=populationSorted(1:20,:);
    
    for n=21:pop_size
        population(n,:)=new_gen(n,:);
    end
  
    %fitness of new population in order to loop again
    for ev=1:pop_size
        population(ev,101)=calc_fit(population(ev,1:100),xy);
    end  
end

%find fittest
fittest=1;
for f=2:pop_size
    if(population(f:101)<population(fittest,101))
        fittest=f;
    end
end

population(fittest,101)=population(fittest,1);

 figure('Name','TSP_GA | Results','Numbertitle','off');
 subplot(2,2,1);
 pclr = ~get(0,'DefaultAxesColor');
 plot(xy(:,1),xy(:,2),'.','Color',pclr);
 title('City Locations');
 subplot(2,2,2);
 rte = population(fittest,1:101);
 plot(xy(rte,1),xy(rte,2),'r.-');
 title(sprintf('Total Distance = %1.4f',getDist(population(fittest,:),xy)));

%Distance Function
function distance_sum=getDist(chrom,xy)

        distance_sum=0;
        for x=2:100
            
            Ax=xy(chrom(x-1),1);
            Ay=xy(chrom(x-1),2);
            
            Bx=xy(chrom(x),1);
            By=xy(chrom(x),2);
            
            d=sqrt( (Bx-Ax)^2 + (By-Ay)^2 );
            
            distance_sum=distance_sum+d;
        end
        %also calculate and add distance from finish city to start city!
        Ax1=xy(chrom(100),1);
        Ay1=xy(chrom(100),2);
            
        Bx1=xy(chrom(1),1);
        By1=xy(chrom(1),2);
        
        d1=sqrt( (Bx1-Ax1)^2 + (By1-Ay1)^2 );
        
        distance_sum=distance_sum+d1;

end

%Fitness function
%Calculate fitness of each chromosome in the population and return an array
%with fitness values
function distance_sum = calc_fit(chrom, xy)


        distance_sum=0;
        for x=2:100
            
            Ax=xy(chrom(x-1),1);
            Ay=xy(chrom(x-1),2);
            
            Bx=xy(chrom(x),1);
            By=xy(chrom(x),2);
            
            d=sqrt( (Bx-Ax)^2 + (By-Ay)^2 );
            
            distance_sum=distance_sum+d;
        end
        %also calculate and add distance from finish city to start city!
        Ax1=xy(chrom(100),1);
        Ay1=xy(chrom(100),2);
            
        Bx1=xy(chrom(1),1);
        By1=xy(chrom(1),2);
        
        d1=sqrt( (Bx1-Ax1)^2 + (By1-Ay1)^2 );
        
        distance_sum=distance_sum+d1;

end

%Turnament selection Function
% fit = array containing the fitness evaluation of each chromosome in pop
function turnament = select_turnmnt(pop_size, fit)
    
    %array to hold the 2 final selected chromosomes
    temp_results=[];
    
%     fitness=sortrows(fit);
%     fitness1=fitness(1:30);
    
    %run 2 turnaments
    for i=1:2
        
        %generate 2 random integers between 1-pop_size
        rands=randi([1,100],2,1);

        %retrieve the fitness evaluation of the 2 chromosomes
        chrom1Fit=fit(rands(1));
        chrom2Fit=fit(rands(2));

        %take the fittest chromosome
        if(chrom1Fit<chrom2Fit)
            temp_results=[temp_results, rands(1)];
        else
            temp_results=[temp_results, rands(2)];
        end  
    end
   turnament=temp_results;
end

%Roulette Selection Function
function roulette = select_rulette(pop_size, fit)

    %array to hold the 2 final selected chromosomes
    temp_array=[];
    
    %calculate propabilities
    summ=sum(fit);
    probabilities=[];
    for i=1:pop_size
        probabilities=[probabilities, fit(i)/summ];
    end
    
    %calculate comulatives
    comulatives=[];
    %append the first comulative (first propability as it is)
    comulatives=[comulatives, probabilities(1)];
    %calculate and append the rest of the comulatives 
    for i=2:pop_size
        comulatives=[comulatives, sum(probabilities(1:i))];
    end
    
    %run rulette 2 times
    for i=1:2
        %generate random number between 0 and 1
        r=rand(1);
        for x=1:pop_size
            if ( comulatives(x) > r)
                temp_array=[temp_array, x];
                r=[];
                break
            end
        end
    end
    r=[];
    roulette=temp_array;
end

%Order Crossover function
function OrderCrossover = Order_crossover(parent1, parent2)

    OrderCrossover=zeros(2,100);
    %generate 2 crossover points
    cross_points=randi([ 1, length(parent1)-1],2,1);
    
    %sequence of bits in parents starting from 2nd cross point
    newSeqParent1=[];
    newSeqParent1=[parent1( (cross_points(2)+1) :100), parent1(1:cross_points(1)), parent1(cross_points(1)+1:cross_points(2))];
    
    newSeqParent2=[];
    newSeqParent2=[parent2( (cross_points(2)+1) :100), parent2(1:cross_points(1)), parent2(cross_points(1)+1:cross_points(2))];
   
    
     %initiate childern with zeros
     child1=zeros(100,1);
     child2=zeros(100,1);
     
     %keep the part between cut point 1 and cut point 2 the same as parents
     child1( cross_points(1)+1 : cross_points(2) ) = parent1( cross_points(1)+1 : cross_points(2) );
     child2( cross_points(1)+1 : cross_points(2) ) = parent2( cross_points(1)+1 : cross_points(2) );
     
     
     %check all genes in child 1
     for i=1:100
         %if allele is 0, it needs to be filed with genes from
         %newSeqParent2
         if (child1(i)==0)
             %Start by testing if 1st element is newSeqParent2 can be
             %used to replace the 0 (chromosomes should not have duplicates
             test=1;
             %do this until the 0 is replaces succesfully and validly
             while(child1(i)==0)
                 %check if this element in newSeqParent2 does not already
                 %exist in child1
                 exist=0;
                 for t=1:100
                     if(newSeqParent2(test)==child1(t))
                         exist=1;
                     end
                 end
                 
                 %if it does not exists, use it to replace the 0 in the
                 %child
                 if( exist==0 )
                     %replace the 0 in child
                     child1(i)=newSeqParent2(test);
                     %remove the used element from newSeqParent2
                     newSeqParent2(test)=[]; 
                 else
                  %try with the next element in newSeqParent2
                  test=test+1;
                 end
                 exist=[];
             end
         end
     end
     

   %check all genes in child 1
     for t=1:100
         %if allele is 0, it needs to be filed with genes from
         %newSeqParent2
         if (child2(t)==0)
             %Start by testing if 1st element is newSeqParent2 can be
             %used to replace the 0 (chromosomes should not have duplicates
             test1=1;
             %do this until the 0 is replaces succesfully and validly
             while(child2(t)==0)
                 %check if this element in newSeqParent2 does not already
                 %exist in child1
                 exist1=0;
                 for y=1:100
                     if(newSeqParent1(test1)==child2(y))
                         exist1=1;
                     end
                 end
                 %if it does not exists, use it to replace the 0 in the
                 %child
                 if( exist1==0 )
                     %replace the 0 in child
                     child2(t)=newSeqParent1(test1);
                     %remove the used element from newSeqParent2
                     newSeqParent1(test1)=[]; 
                 else
                  %try with the next element in newSeqParent2
                  test1=test1+1;
                 end
                 exist1=[];
             end
         end
     end
     

    OrderCrossover(1,:)=child1;
    OrderCrossover(2,:)=child2;
    
end

%Partial Mapping crossover Function
function PartMapCrossover = PartMap_crossover(parent1, parent2)

    PartMapCrossover=[];
    %generate 2 crossover points
    cross_points=randi([ 1, length(parent1)-1],2,1);
    
    %mapping system
    mapping1=[];
    mapping1=parent1(cross_points(1)+1 : cross_points(2));
    
    mapping2=[];
    mapping2=parent2(cross_points(1)+1 : cross_points(2));
    

     %initiate childern with zeros
     child1=zeros(100,1);
     child2=zeros(100,1);
     
     %keep the part between cut point 1 and cut point 2 the same as OTHER
     %parent
     child1( cross_points(1)+1 : cross_points(2) ) = parent2( cross_points(1)+1 : cross_points(2) );
     child2( cross_points(1)+1 : cross_points(2) ) = parent1( cross_points(1)+1 : cross_points(2) );
     
     %To fill the rest of the genes for each child, put the genes from
     %parent1 to child1 at the same position (check that value does not
     %exists already). If the value already exists, use the mapping2 for
     %child1 etc.
     
      %check all genes of child
      for k=1:100
          %if its 0, it needs to be filled
          if(child1(k) == 0)
              %check if the parent gene on the same position can be used
              exist=0;
              for l=1:100
                  if (parent1(k) == child1(l))
                        exist=1;
                  end
              end
              if(exist==0)
                  child1(k)=parent1(k);
              else
                  %if not, try to use the mapping of the other parent
                  while(child1(k)==0)
                      %try the first value of the mapping
                      tryit2=1;
                      %check if it exists
                      exist2=0;
                      for m=1:100
                          if(mapping2(tryit2)==child1(m))
                              exist2=1;
                          end
                      end
                      %use it if its valid
                      if(exist2==0)
                          child1(k)=mapping2(tryit2);  
                          %remove it from the mapping
                          mapping2(tryit2)=[];
                      else
                          %go to the next element on the mapping to try it
                          tryit2=tryit2+1;
                      end
                  end
              end
          end
      end
      

    
    for l=1:100
          %if its 0, it needs to be filled
          if(child2(l) == 0)
              %check if the parent gene on the same position can be used
              exist1=0;
              for m=1:100
                  if (parent2(l) == child2(m))
                        exist1=1;
                  end
              end
              if(exist1==0)
                  child2(l)=parent2(l);
              else
                  %if not, try to use the mapping of the other parent
                  while(child2(l)==0)
                      %try the first value of the mapping
                      tryitt=1;
                      %check if it exists
                      existt=0;
                      for n=1:100
                          if(mapping1(tryitt)==child2(n))
                              existt=1;
                          end
                      end
                      %use it if its valid
                      if(existt==0)
                          child2(l)=mapping1(tryitt);  
                          %remove it from the mapping
                          mapping1(tryitt)=[];
                      else
                          %go to the next element on the mapping to try it
                          tryitt=tryitt+1;
                      end
                  end
              end
          end
    end
      

      PartMapCrossover=[child1, child2];
    
end

%Swap Mutation Function
function chromosomes = swapMutation(chrom1, chrom2)

    chromosomes=zeros(2,100);

    %generate 2 random numbers between 1 and 100
    swap_points=randi([ 1, 100],2,1);
    
    %do the swapping
    temp1=chrom1(swap_points(1));
    temp2=chrom1(swap_points(2));
    
    chrom1(swap_points(1))=temp2;
    chrom1(swap_points(2))=temp1;
    
    %do the swapping
    temp1=chrom2(swap_points(1));
    temp2=chrom2(swap_points(2));
    
    chrom2(swap_points(1))=temp2;
    chrom2(swap_points(2))=temp1;
    
    chromosomes(1,1:100)=chrom1;
    chromosomes(2,1:100)=chrom2;
    
   
end

%Flip Mutation Function
function chromosome = flipMutation(chrom)

    %generate 2 random numbers between 1 and 100
    swap_points=randi([ 1, 100],2,1);
    
    a=swap_points(1);
    b=swap_points(2);
    
    %do the flipping
    cont=1;
    while(cont==1)
        temp1=chrom(a);
        temp2=chrom(b);
        
        chrom(b)=temp1;
        chrom(a)=temp2;
        
        if(a<b)
            a=a+1;
            b=b-1;
        else
            a=a-1;
            b=b+1;
        end
        
        if(a==b)
            cont=0;
        end
            
    end    
    
    chromosome=chrom;
    
end








    