function [B,V_A,U_A] = Bidiag_Francis_Step_Update_U_V(B,V_A,U_A)
    [m, n] = size(B);

   if m==2
            T11 = B(1,1);
            T21 = B(1,2);
            Tmm = B(m,m);

         G = Givens_rotation( [(T11-Tmm)
                                         T21      ]);
        
         B=B*G;
        display(B);
        B(1,2)=0;
        B(2,1)=0;
       

            T211=V_A(1,1)^2;
             T221=V_A(1,2)*V_A(1,1);
             T2mm=V_A(m-1,m)^2 + V_A(m,m)^2;

             T311=U_A(1,1)^2;
             T321=V_A(1,2)*V_A(1,1);
             T3mm=V_A(m,m);

            G2 = Givens_rotation( [(T211-T2mm)
                                         T221      ]);
            G3 = Givens_rotation( [(T311-T3mm)
                                         T321      ]);
         V_A=V_A*G2;
        display(V_A);
        V_A(1,2)=0;
        V_A(2,1)=0;

        U_A=U_A*G3;
        display(U_A);
        U_A(1,2)=0;
        V_A(2,1)=0;
   return
   end



  
    for i = 1:min(m-1, n-1)
            
  
        %last statement
        if i == m-1
             T11=B(m-1,m-1);
             T21=B(m-1,m);
             Tmm=B(m,m);
              G = Givens_rotation( [(T11-Tmm)
                                         T21      ]);
             
             % Apply the rotation to the (m-1,m) and (m,m) elements
             B(m-1:m, m-1:m) = G' *B(m-1:m, m-1:m);
             B(m,m-1) = 0;

             T211=V_A(m-1,m-1);
             T221=V_A(m-1,m);
             T2mm=V_A(m,m);
             G2 = Givens_rotation( [(T211-T2mm)
                                         T221      ]);
            
             V_A(m-1:m, m-1:m)= G2'*V_A(m-1:m, m-1:m);

             T311=U_A(m-1,m-1);
             T321=V_A(m-1,m);
             T3mm=V_A(m,m);

              G3 = Givens_rotation( [(T311-T3mm)
                                         T321      ]);
             U_A(m-1:m, m-1:m)= G3'*U_A(m-1:m, m-1:m);

             display(B);
             display(U_A);
             display(V_A);
            
            
           
             return
        end
 
        if(i == 1)
            %compute G0 and add bulge 
            T11 = B(1,1)^2;
            T21 = B(1,2)*B(1,1);
            Tmm = B(m-1,m)^2 + B(m,m)^2;

             T211=V_A(1,1)^2;
             T221=V_A(1,2)*V_A(1,1);
             T2mm=V_A(m-1,m)^2 + V_A(m,m)^2;

             T311=U_A(1,1)^2;
             T321=V_A(1,2)*V_A(1,1);
             T3mm=V_A(m,m);
            
           

            % Compute the Givens rotation that annihilates the (1,2) and (2,1) elements
            G = Givens_rotation( [(T11-Tmm)
                                         T21      ]);
            G2 = Givens_rotation( [(T211-T2mm)
                                         T221      ]);
            G3 = Givens_rotation( [(T311-T3mm)
                                         T321      ]);
            %display(G);
            % Apply the rotation to the first two 2X2 matrix
            B(1:2, 1:2) = B(1:2, 1:2) * G;
            U_A(1:2,1:2) =U_A(1:2,1:2) *G2;
            V_A(1:2,1:2) =V_A(1:2,1:2) *G3;
            display(B);
        end

        % removed buldge on bottom
       
        T11 = B(i, i);
        T21 = B(i+1,i);
        Tmm = B(m,m);

         T211=V_A(i,i);
         T221=V_A(i+1,i);
         T2mm=V_A(m,m);

         T311=U_A(i,i);
         T321=U_A(i+1,i);
         T3mm=U_A(m,m);
        
          G = Givens_rotation( [T11- Tmm
                                     T21      ]);
         G2 = Givens_rotation( [T211- T2mm
                                     T221      ]);
         G3 = Givens_rotation( [T311- T3mm
                                     T321      ]);
        %add bulge and set the previous bulge to 0 
        B(i:i+1, i:i+2) = G' * B(i:i+1,i:i+2);
        B(i+1,i) = 0;
        display(B);

        V_A(i:i+1, i:i+2) = G2' * V_A(i:i+1,i:i+2);
        V_A(i+1,i) = 0;
        display(V_A);
        
        U_A(i:i+1, i:i+2) = G3' * U_A(i:i+1,i:i+2);
        U_A(i+1,i) = 0;
        display(U_A);

        % removed buldge 
        T11 = B(i, i+1);
        T21 = B(i, i+2);
        Tmm = B(m,m);

        T211 = V_A(i, i+1);
        T221 = V_A(i, i+2);
        T2mm = V_A(m,m);

        T311 = U_A(i, i+1);
        T321 = U_A(i, i+2);
        T3mm = U_A(m,m);

        

      G = Givens_rotation( [T11- Tmm
                                     T21      ]);
         G2 = Givens_rotation( [T211- T2mm
                                     T221      ]);
         G3 = Givens_rotation( [T311- T3mm
                                     T321      ]);
        %display(G);
        B(i:i+2,i+1:i+2) =  B(i:i+2,i+1:i+2)*G ;
        display(B);
        B(i,i+2) = 0;
        display(B);
        B(m-1,m)=0;

         V_A(i:i+2,i+1:i+2) =  V_A(i:i+2,i+1:i+2)*G2 ;
        display(B);
        V_A(i,i+2) = 0;
        display(V_A);
        V_A(m-1,m)=0;

           U_A(i:i+2,i+1:i+2) =  U_A(i:i+2,i+1:i+2)*G3 ;
        display(B);
        U_A(i,i+2) = 0;
        display(U_A);
        U_A(m-1,m)=0;

    end

end