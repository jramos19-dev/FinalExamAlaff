function B = Bidiag_Francis_Step(B)
    [m, n] = size(B);

   if m==2
            T11 = B(1,1)^2;
            T21 = B(1,2)*B(1,1);
            Tmm = B(m-1,m)^2 + B(m,m)^2;

         G = Givens_rotation( [(T11-Tmm)
                                         T21      ]);
        
         B=B*G;
        display(B);
        B(1,2)=0;
        B(2,1)=0;
       
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
             display(B);
            
            
           
             return
        end
 
        if(i == 1)
            %compute G0 and add bulge 
            T11 = B(1,1)^2;
            T21 = B(1,2)*B(1,1);
            Tmm = B(m-1,m)^2 + B(m,m)^2;

            % Compute the Givens rotation that annihilates the (1,2) and (2,1) elements
            G = Givens_rotation( [(T11-Tmm)
                                         T21      ]);
            %display(G);
            % Apply the rotation to the first two 2X2 matrix
            B(1:2, 1:2) = B(1:2, 1:2) * G ;
            
            display(B);
        end

        % removed buldge on bottom
       
        T11 = B(i, i);
        T21 = B(i+1,i);
        Tmm = B(m,m);
        
        
        G = Givens_rotation( [T11- Tmm
                                     T21      ]);

        %add bulge and set the previous bulge to 0 
        B(i:i+1, i:i+2) = G' * B(i:i+1,i:i+2);
        B(i+1,i) = 0;
        display(B);

        % removed buldge 
        T11 = B(i, i+1);
        T21 = B(i, i+2);
        Tmm = B(m,m);

        % calculate the givens for top right 
        G = Givens_rotation( [T11- Tmm
                                     T21      ]);
        %display(G);
        B(i:i+2,i+1:i+2) =  B(i:i+2,i+1:i+2)*G ;
        display(B);
        B(i,i+2) = 0;
        display(B);
        B(m-1,m)=0;


    end

end