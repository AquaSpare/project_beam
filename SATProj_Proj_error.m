


[u,v] = beam_ana();

iter = 8;
N = 15;
errorProj = zeros(1,iter);
errorSAT = zeros(1,iter);
stepsProj = zeros(1,iter);
stepsSAT = zeros(1,iter);


for i = 1:iter
<<<<<<< HEAD
    
    [solProj,tProj,xProj,hProj,kProj] = beam_eq_projection((N*i)+1,-1,0,1,0.14,1,4,4,2,100,0.0001,{u,v});
=======
    [solProj,tProj,xProj,hProj,kProj] = beam_eq_projection_dav((N*i)+1,-1,0,1,0.14,1,4,4,2,0.0001,{u,v});
>>>>>>> 5e98781d4356a2722b6a702013f02998bfd5a830
    exact = [u(xProj(1:(N*i)+1),tProj(end)) v(xProj((N*i)+2:end),tProj(end))];
    stepsProj(i) = hProj;
    errorProj(i) = sqrt(1/(N*i))*norm(exact'-solProj(:,end),2);
                                        
    [solSAT, tSAT, xSAT, hSAT, kSAT] = SAT_projection_dav((N*i)+1,-1,0,1,0.14, 1, 4, 4, 2, 100,0.0001,{u,v});
    exact = [u(xSAT(1:(N*i)+1),tSAT(end)) v(xSAT((N*i)+2:end),tSAT(end))];
    stepsSAT(i) = hSAT;
    errorSAT(i) = sqrt(1/(N*i))*norm(exact'-solSAT(:,end),2);
    
end
figure(1)
loglog(stepsProj,errorProj,'r*')
hold on
loglog(stepsSAT,errorSAT,'b*')
xlabel('log(h)')
ylabel('L2 norm of error')
legend('Projection', 'SAT');