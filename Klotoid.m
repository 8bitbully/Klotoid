classdef Klotoid
    properties (Dependent = true, Access = public)
        TegetinAcisi % Ks noktasindaki tegetin acisi. To(grad)
        Parametre % A
        DaireninYaricapi % R
        KlotoidinBoyu % L
        dikKoordinatlar % Ks'nin dik koordinatlari [X, Y]
        daireMerkeziKoordinatlari % [Xm, Ym]
        RakordmanPayi % dR
        Tegetler % [Tk, Tu] : Tk-kisa teget, Tu-uzun teget.
        kutupsalKoordinatlar % Ks' nin kutupsal koordinatlari [S, sigma] : sigma(grad)
    end
    properties (Access = private)
        Delta % delta acisi
        DaireninMerkezi % M
        A
        R
        L
        T
        dikKoor
        kutupsalKoor
        daireMerkezi
        teget
        deltaR
        D % delta
    end
 
    methods
        function this = Klotoid(R, L, A)
            format longG
            this.DaireninYaricapi = R;
            this.KlotoidinBoyu = L;
            this.Parametre = A;
        end
    end
    
    % GET METHODS
    methods
        function a = get.Parametre(this)
            a = this.A;
        end
        function r = get.DaireninYaricapi(this)
            r = this.R;
        end
        function l = get.KlotoidinBoyu(this)
            l = this.L;
        end
        function t = get.TegetinAcisi(this)
            t = (this.T);
        end
        function rp = get.RakordmanPayi(this)
            rp = this.deltaR;
        end
        function x = get.dikKoordinatlar(this)
            x = this.dikKoor;
        end
        function kc = get.kutupsalKoordinatlar(this)
            kc = this.kutupsalKoor;
        end
        function d = get.Delta(this)
            d = this.D;
        end
        function xm = get.daireMerkeziKoordinatlari(this)
            xm = this.daireMerkezi;
        end
        function tg = get.Tegetler(this)
            tg = this.teget;
        end
    end
    
    % SET METHODS
    methods
        function this = set.DaireninYaricapi(this, r)
            this.R = r;
        end
        function this = set.KlotoidinBoyu(this, l)
            this.L = l;
            r = this.R;
            t = l * (200 / pi) / (2*r);
            this.T = round(t, 5);
        end
        function this = set.Parametre(this, a)
            l = this.L;
            r = this.R;
            this.A = a;
            t = this.T*pi/200;
            
            x = l - l^5 / (40*a^4);
            y = l^3 / (6*a^2) - l^7 / (336*a^6);
            Ym = y + r * cos(t);
            Xm = x - r * sin(t);
            dR = Ym - r;
            Tk = y / sin(t);
            Tu = x - y * cot(t);
            s = sqrt(x^2+y^2);
            sigma = (atan(y / x))*200/pi;
            
            this.deltaR = round(dR, 5);
            this.dikKoor = round([x, y], 3);
            this.daireMerkezi = round([Xm, Ym], 3);
            this.teget = round([Tk, Tu], 3);
            this.kutupsalKoor = [round(s, 3), round(sigma, 5)];
        end
    end
    
    methods
        function [t_, T_, alpha, b_, BS] = deltaTegetUzunlugu(klotoid, delta)
            R_ = klotoid.R;
            dR_ = klotoid.deltaR;
            T_ = klotoid.T;
            Xm = klotoid.daireMerkeziKoordinatlari(1);
            
            t_ = round((R_ + dR_)*tan((delta*pi/200) / 2), 3); % delta: grad to radians
            alpha = round(delta - 2*T_, 5);
            b_ = round(R_*pi*alpha/200, 3);
            MS = t_ / sin((delta*pi/200) / 2);
            BS = round(MS - R_, 3);
            T_ = round(t_ + Xm, 3);
        end
        
        function [T_, Z_, BS_] = kurbunTamamiKlotoid(klotoid)
            x = klotoid.dikKoordinatlar(1);
            y = klotoid.dikKoordinatlar(2);
            To_ = klotoid.T;
            
            Z_ = round(y*tan(To_ * pi / 200),3);
            T_ = round(Z_ + x, 3);
            BS_ = round(y / cos(To_ * pi / 200), 3);
        end
        
        % KB1 = QS - T
        function [KS1, B, KS2, KB2, app] = kilometraj(klotoid, KB1, b, aralik)
            l = klotoid.L;
            
            KS1 = KB1 + l;
            B = KS1 + b/2;
            KS2 = KB1 + l + b;
            KB2 = KB1 + 2*l + b;
            
            counter = 1;
            points = [];
            for i = KB1-mod(KB1,aralik)+aralik:aralik:KB2
                points.pointer(counter) = i;
                counter = counter + 1;
            end
            point = sort([points.pointer, KB1, KS1, B, KS2, KB2]);
            
            for i = 1:length(point)
                r = Klotoid.findR(point(i)-point(1), klotoid.A);
                params.App(i) = Klotoid(r, point(i)-point(1), klotoid.A);
            end
            
            app.km = point;
            app.Params = [params.App];
        end
    end
    
    
    %Static Method
    methods (Static)
        function A = findA(R, L)
            A = round(sqrt(R * L), 3);
        end
        
        function R = findR(L, A)
            R = round(A * A / L, 3);
        end
        
        function L = findL(R, A)
            L = round(A *A / R, 3);
        end
        % To_: grad, L_: metre
        function [To_, L_] = delta2L(delta, R)
            To_ = delta / 2;
            L_ = round((To_ * 2 * R) / (200/pi), 3);
        end
        
        % To_: grad, R_ metre
        function [To_, R_] = delta2R(delta, L)
            To_ = delta / 2;
            R_ = round(( L*(200/pi) ) / (2*To_), 3);
        end
        
        function params = findMainParams(app, aralik)
            index = mod(app.km, aralik) ~= 0;
            params.Parameter = app.Params(index);
            params.km = app.km(index);
        end
    end
end