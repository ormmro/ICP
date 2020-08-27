laser_A = load('laser_A.txt');
laser_B = load('laser_B.txt');
A_hat = zeros(size(laser_A));
d_B = zeros(size(laser_B));
w = zeros(2,2);
n = 0;
E = [0 0];
while n<20
    if n > 0
        laser_A = A_hat;
    end
    
    x_A = laser_A(:,1);
    y_A = laser_A(:,2);
    x_B = laser_B(:,1);
    y_B = laser_B(:,2);
    
    %find new laser_B data
    for a=1:360
        for b = 1:360
            dist = sqrt(power((x_A(a)-x_B(b)),2) +power((y_A(a)-y_B(b)),2));
            if b == 1 || dist < distOld
                num_short = b;
                distOld = dist;
            end 
        end
        d_B(a,:) = laser_B(num_short,:);
    end 

    %find center of mass
    center_A(:,1) = sum(x_A)/length(x_A);
    center_A(:,2) = sum(y_A)/length(y_A);
    center_d_B(:,1) = sum(d_B(:,1))/length(d_B(:,1));
    center_d_B(:,2) = sum(d_B(:,2))/length(d_B(:,2));

    %set every point to same center
    d_A = laser_A - center_A;
    dd_B = d_B - center_d_B;
    
    %find rotation and translation
    for i = 1:360
        w = w + d_A(i,:)'*dd_B(i,:);
    end
    [U,S,V] = svd(w);
    R = V*U';
    T = (center_A' - (R*center_d_B'));

    %change rotation and translation    
    for i = 1:360
        laser_A_tran = [laser_A(i,:)'; 1];
        hat = [R -T;0 0 1]*laser_A_tran;
        A_hat(i,:) = hat(1:2,:)';
    end
    
    %find error
    for i = 1:360
        E = E + abs((A_hat(i,:) - d_B(i,:)));
    end
    n = n + 1;
    w=0;
end
%plot
figure(1)
plot(x_B,y_B,'.','color','red')
hold on
plot(A_hat(:,1),A_hat(:,2),'.','color','green')

