classdef AdaBoost
    methods
        function F = CreateFeatures(obj, input, sketch_count, end_ind_list)
%             x = [0 end_ind_list];
% %             F = obj.GetRubineFeatures(input);
%             f1 = sketch_count;
%             [f2 bb] = obj.GetAspectRatio(input);
%             f3 = obj.GetPHist(bb, input);
%             [f4, angles] = obj.GetAHist(input(:,1:2));
%             f5 = obj.GetFLDist(sketch_count,input,x);
%             f6 = obj.GetArcLength(sketch_count,input,x);
%             f7 = obj.GetStrokeLineBool(sketch_count,input(:,1:2),x);
%             f8 = 0; %obj.GetMaxAngle(sketch_count,angles,x);
%             f9 = obj.GetStrokeArea(sketch_count,input(:,1:2),x);
%             f10 = obj.GetSideRatios(sketch_count, input, x, bb);
%             f11 = obj.GetTBRatios(sketch_count, input, x, bb);
%             f12 = 0;% obj.GetMinMaxFeatures(sketch_count, input, x);
%             F = [f1;f2;f3;f4;f5;f6;f9;f10;f11];
%             s_f1 = size(f1);
%             s_f2 = size(f2);
%             s_f3 = size(f3);
%             s_f4 = size(f4);
%             s_f5 = size(f5);
%             s_f6 = size(f6);
%             s_f7 = size(f7);
%             s_f8 = size(f8);
%             s_f9 = size(f9);
%             s_f10 = size(f10);
%             s_f11 = size(f11);
%             s_f12 = size(f12);
            F = obj.GetRubineFeatures(input);
            s_F = size(F)
        end
        
        function [f2, bb] = GetAspectRatio(obj, input)
           min_x =  min(input(:,1));
           min_y =  min(input(:,2));
           max_x =  max(input(:,1));
           max_y =  max(input(:,2));
           bb = [min_x min_y; max_x max_y];
           f2 = (max_x-min_x)/(max_y-min_y);
        end
        
        function f3 = GetPHist(obj, bbox, input)
            len = bbox(2,2) - bbox(1,2);
            wdt = bbox(2,1) - bbox(1,1);
            grid_len = len/3;
            grid_wdt = wdt/3;
            gl = [bbox(1,2); bbox(1,2) + grid_len; bbox(1,2) + grid_len*2; bbox(1,2) + grid_len*3];
            gw = [bbox(1,1); bbox(1,1) + grid_wdt; bbox(1,1) + grid_wdt*2; bbox(1,1) + grid_wdt*3];
            f3 = zeros(3);
            for i = 1:3
                for j = 1:3
                    a = find((input(:,1)>=gw(i) & input(:,1)<gw(i+1)) & (input(:,2)>=gl(j) & input(:,2)<gl(j+1)));
                    f3(i,j) = length(a);
                end
            end
            for i = 1:3
                a = find((input(:,1)>=gw(i) & input(:,1)<gw(i+1)) & (input(:,2)==gl(4)));
                f3(i,3) = f3(i,3) + length(a);
                b = find((input(:,1)==gw(4)) & (input(:,2)>=gl(i) & input(:,2)<gl(i+1)));
                f3(3,i) = f3(3,i) + length(b);
            end
            f3 = reshape(f3,[9,1]);
        end
        
        function [f4, angles] = GetAHist(obj,input)
            d = diff(input);
            theta = atan2d(d(:,2), d(:,1));
            f4 = zeros(8,1);
            for i = 1:4
               f4(i) = length(find(theta >= (i-1)*45 & theta < i*45));
               f4(i+4) = length(find(theta < (i-1)*(-45) & theta >= i*(-45)));
            end
            f4(8) = f4(8) + length(find(theta == 180));
            angles = theta;
        end
        
        function f5 = GetFLDist(obj, sketch_count, input, end_ind_list)
            d = diff(input);
            %dist = sqrt(d(:,1).^2 + d(:,2).^2);
            f5 = zeros(sketch_count,1);
            for i = 1:sketch_count
                dist = input(end_ind_list(i+1),:) - input(end_ind_list(i) + 1,:);
                f5(i) = sqrt(dist(1)^2 + dist(2)^2);
                %f5(i) = dist(end_ind_list(i) + 1);
            end
            f5 = mean(f5);
        end
        
        function f6 = GetArcLength(obj, sketch_count, input, end_ind_list)
            d = diff(input);
            dist = sqrt(d(:,1).^2 + d(:,2).^2);
            f6 = 0;
            for i = 1:sketch_count
                indf = end_ind_list(i) + 1;
                indl = end_ind_list(i+1) - 1;
                f6 = f6 + sum(dist(indf:indl));
            end
        end
        
        function f7 = GetStrokeLineBool(obj, sketch_count, input, end_ind_list)
            f7 = zeros(sketch_count,1);
            for i = 1:sketch_count
                p1 = input(end_ind_list(i)+1, :);
                p2 = input(end_ind_list(i+1), :);
                p = input(end_ind_list(i)+1:end_ind_list(i+1),:);
                a = p1 - p2;
                b = p - p2;
                z = zeros(length(b),1);
                b = [b z];
                c = ones(length(b),2);
                c = c .* a;
                c = [c z];
                d = abs(cross(c,b));
                max_d = max(d(:,3));
                e = find(d(:,3) >= 0.0150);
                if(length(e) >= 5)
                    f7(i) = 0;
                else
                    f7(i) = 1;
                end
            end
        end
        
        function f8 = GetMaxAngle(obj, sketch_count, input, end_ind_list)
            f8 = zeros(sketch_count,1);
            for i = 1:sketch_count
                f8(i) = max(input(end_ind_list(i)+1:end_ind_list(i+1)-1));
            end
        end
        
        function f9 = GetStrokeArea(obj, sketch_count, input, end_ind_list)
            v = input - input(1,:);
            f9 = zeros(sketch_count,1);
            for i = 1:sketch_count
                startp = end_ind_list(i) + 1;
                lastp = end_ind_list(i+1);
                v = input(startp+1:lastp, :) - input(startp, :);
                v1 = v(1:end-1, :);
                v2 = v(2:end, :);
                z = zeros(length(v1),1);
                v1 = [v1 z];
                v2 = [v2 z];
                vcp = cross(v1,v2);
                vcp = vcp(:,3)/2;
                sgnvcp = sign(vcp);
                s = vcp.*sgnvcp;
                f9(i) = sum(s);
            end
            avg_s = mean(f9);
            f9 = avg_s;
        end
        
        function f10 = GetSideRatios(obj, sketch_count, input, end_ind_list, bb)
            f10 = zeros(sketch_count,2);
            for i = 1:sketch_count
                startp = input(end_ind_list(i)+1,1);
                lastp = input(end_ind_list(i+1),1);
                f10(i,1) = (startp - bb(1,1))/(bb(2,1) - bb(1,1));
                f10(i,2) = (lastp - bb(1,1))/(bb(2,1) - bb(1,1));
            end
            f101 = mean(f10(:,1));
            f102 = mean(f10(:,2));
            f10 = [f101; f102];
        end
        
        function f11 = GetTBRatios(obj, sketch_count, input, end_ind_list, bb)
            f11 = zeros(sketch_count,2);
            for i = 1:sketch_count
                startp = input(end_ind_list(i)+1,2);
                lastp = input(end_ind_list(i+1),2);
                f11(i,1) = (bb(2,2) - startp)/(bb(2,2) - bb(1,2));
                f11(i,2) = (bb(2,2) - lastp)/(bb(2,2) - bb(1,2));
            end
            f111 = mean(f11(:,1));
            f112 = mean(f11(:,2));
            f11 = [f111; f112];
        end
        
        function f12 = GetMinMaxFeatures(obj,sketch_count, input, end_ind_list)
            f12 = zeros(sketch_count,10);
            for i =1:sketch_count
                startp = end_ind_list(i) + 1;
                lastp = end_ind_list(i+1);
                d = diff(input(startp:lastp,:));
                x = input(startp:lastp,1);
                y = input(startp:lastp,2);
                x_dir = d(:,1);
                y_dir = d(:,2);
                x_dir = sign(x_dir);
                y_dir = sign(y_dir);
                
                x_dir_diff = diff(x_dir);
                y_dir_diff = diff(y_dir);
                
                id_x = find(abs(x_dir_diff) >=1);
                id_y = find(abs(y_dir_diff) >=1);
                
                lmin_x = length(find(x_dir_diff==2));
                lmax_x = length(find(x_dir_diff==-2));
                x_start_dir = x_dir(1);
                x_last_dir = x_dir(end);
                x_last = x(id_x(end));
                x_dist = x(end) - x_last;
                
                
                lmin_y = length(find(y_dir_diff==2));
                lmax_y = length(find(y_dir_diff==-2));
                y_start_dir = y_dir(1);
                y_last_dir = y_dir(end);
                y_last = y(id_y(end));
                y_dist = y(end) - y_last;
                f12(i,:) = [lmin_x, lmax_x, x_start_dir, x_last_dir, x_dist, lmin_y, lmax_y, y_start_dir, y_last_dir, y_dist];
            end
            f121 = sum(f12(:,1));
            f122 = sum(f12(:,2));
            a = length(find(f12(:,3) == -1));
            b = length(find(f12(:,4) == -1));
            f123 = 0;
            f124 = 0;
            if(a >=5)
               f123 = -1;
            else
                f123 = 1;
            end
            if(b >=5)
               f124 = -1;
            else
                f124 = 1;
            end
            f125 = mean(f12(:,5));
            
            f126 = sum(f12(:,6));
            f127 = sum(f12(:,7));
            a = length(find(f12(:,8) == -1));
            b = length(find(f12(:,9) == -1));
            f128 = 0;
            f129 = 0;
            if(a >=5)
               f128 = -1;
            else
               f128 = 1;
            end
            if(b >=5)
               f129 = -1;
            else
               f129 = 1;
            end
            f1210 = mean(f12(:,10));
            f12 = [f121;f122;f123;f124;f125;f126;f127;f128;f129;f1210];
        end
        
        function F = GetRubineFeatures(obj, input)
            points = input;
            x31 = points(3,1) - points(1,1);
            y31 = points(3,2) - points(1,2);
            f1 = x31/sqrt((x31^2) + (y31^2));
            f2 = y31/sqrt((x31^2) + (y31^2));

            points_max = max(points);
            points_min = min(points);
            BB = [points_min(1,1) points_min(1,2);points_max(1,1) points_max(1,2)];
            points_maxmin = points_max - points_min;
            f3 = (points_maxmin(1,1)^2 + points_maxmin(1,2)^2)^(1/2);
            f4 = atan(points_maxmin(1,2)/points_maxmin(1,1));

            xend1 = points(end,1) - points(1,1);
            yend1 = points(end,2) - points(1,2);
            f5 = (xend1^2 + yend1^2)^(1/2);
            f6 = xend1/f5;
            f7 = yend1/f5;

            dx = (points(2:end,1) - points(1:end - 1,1));
            dy = (points(2:end,2) - points(1:end - 1,2));
            d = sqrt((dx.^2)+(dy.^2));
            f8 = sum(d);

            theta_ptop = (dx(2:end).*dy(1:end-1)) - (dx(1:end-1).*dy(2:end));
            theta_pbottom = (dx(2:end).*dx(1:end-1)) - (dy(2:end).*dy(1:end-1));
            theta_p = atand(theta_ptop/theta_pbottom);
            f9 = sum(theta_p,'all');
            f10 = sum(abs(theta_p),'all');
            f11 = sum(theta_p.^2,'all');

            dt = (points(2:end,3) - points(1:end - 1,3));
            s = (dx.^2 + dy.^2)./(dt.^2);
            f12 = max(s);
            f13 = sum(points(:,3));

            %{
            sizef1 = size(f1);
            sizef2 = size(f2);
            sizef3 = size(f3);
            sizef4 = size(f4);
            sizef5 = size(f5);
            sizef6 = size(f6);
            sizef7 = size(f7);
            sizef8 = size(f8);
            sizef9 = size(f9);
            sizef10 = size(f10);
            sizef11 = size(f11);
            sizef12 = size(f12);
            sizef13 = size(f13);
            %}
            F = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11;f12;f13];
        end
        
        function [alphas, means, pairs] = TrainModel(obj, F)
            y1 = ones(10,1);
            y2 = (ones(10,1)).*(-1);
            y = [y1;y2]';
            alphas = [];
            means = [];
            pairs = [];
            for i = 0:9
                f1 = F(:,:,i+1);
                for j = 0:9
                    if(i<j)
                        pair = [i,j];
                        alpha = [];
                        mean = [];
                        f2 = F(:,:,j+1);
                        H = 0;
                        dt = ones(1,20)./20;
                        dt1 = zeros(1,20);
                        for t = 1:(1000*13)
                            k = mod(t-1,13) + 1;
                            u1 = (f1(k,:) * dt(1:10)');
                            u2 = (f2(k,:) * dt(11:20)');
                            f = [f1(k,:) f2(k,:)];
                            m1 = abs(f - u1);
                            m2 = abs(f - u2);
                            h = (m1 <= m2);
                            h(h==0) = -1;
                            b = (h~=y);
                            if(sum(b) > 10)
                                h(h==-1) = 0;
                                h = imcomplement(h);
                                h(h==0) = -1;
                                b = (h~=y);
                            end
                            e = sum(b .* dt);
                            a = log((1-e)/e)/2;
                            dt1 = dt.*exp(-((a.*y).*h));
                            dt = dt1./(sum(dt1));
                            H = H + (a.*h);
                            alpha = [alpha; a];
                            mean = [mean ;[u1,u2]];
                        end
                        H = sign(H);
                        pairs = [pairs; pair];
                        alphas = [alphas alpha];
                        means = cat(3,means,mean);
%                         for k = 1:30
%                             u1 = zeros(13,1);
%                             u2 = zeros(13,1);
%                             a = zeros(13,1);
%                             e = zeros(13,1);
%                             hT = 0;
%                             for m = 1:13
%                                 u1(m) = (f1(m,:)*d(1:10)')/(sum(d(m,1:10)));
%                                 u2(m) = (f2(m,:)*d(11:20)')/(sum(d(m,11:20)));
%                                 um = (u1(m)+u2(m))/2;
%                                 f = [f1(m,:) f2(m,:)];
%                                 h = d(m,:) .* f;
%                                 b = (h >= um);
%                                 b(b==0) = -1;
%                                 c = (b~=y);
%                                 if(sum(c)>10)
%                                     c = (c-1).*(-1);
%                                 end
%                                 e(m) = sum(c.*d(m,:));
%                                 a(m) = (log((1-e(m))/e(m)))/2;
%                                 z = 1;
%                                 if(m ~= 13)
%                                     size(y);
%                                     size(h);
%                                     d(m+1,:) = d(m,:).*(exp((a(m).*y).*h));
%                                     z = sum(d(m+1,:));
%                                     d(m+1,:) = d(m+1,:)./z;
%                                 else
%                                     d(1,:) = d(m,:).*(exp((a(m).*y).*h));
%                                     z = sum(d(1,:));
%                                     d(1,:) = d(1,:)./z;
%                                 end
%                                 hT = hT + (a(m).*h);
%                             end
%                             H = H + hT;
%                             a = a;
%                             e = e;
%                         end
                        %H = sign(H);
                    end
                end
            end
        end
        
        function class = Classify(obj, F, alphas, means, pairs)
            x = [0,1,2,3,4,5,6,7,8,9];
            b = size(alphas);
            sgnHs = [];
            Hs = [];
            for i = 1:b(2)
                H = 0;
                alpha = alphas(:,i);
                mean = means(:,:,i);
                for j = 1:(b(1))
                    k = mod(j-1,13) + 1;
                    a = alpha(j);
                    u1 = mean(j,1);
                    u2 = mean(j,2);
                    m1 = abs(F(k) - u1);
                    m2 = abs(F(k) - u2);
                    h = (m1 <= m2);
                    H = H + (a*h);
                end
                Hs = [Hs H];
            end
            Hs = Hs;
            sgnHs = sign(Hs);
%             pairs = [2,3;2,7;3,7];
            G = sgnHs;
            G(sgnHs==-1) = 2;
            G(sgnHs==0) = 1;
            inds = zeros(length(G),1);
            for m = 1:length(G)
               inds(m) =  pairs(m,G(m));
            end
            inds = inds;
            count = [];
            
            for n = 0:9
                c = find(inds == n);
                count = [count length(c)];
            end
            c_max = max(count);
            abc = find(count == c_max);
            abc = abc - 1;
            if(length(abc) > 1)
                d = [];
                indexes = [];
                afg = 0;
                for p = 1:length(abc)
                    c = find(inds == abc(p));
                    for q = 1:length(c)
                        d = [d Hs];
                    end
                    afg = afg + length(c);
                    indexes = [indexes afg];
                end
                d_max = max(d);
                d_m_ind = find(d == d_max);
                hm = find(indexes <= d_m_ind);
                hm(end)
                class = abc(hm(end));
            else
                class = x(find(count == c_max));
            end
%             class = class + 1;
        end
        
    end
end