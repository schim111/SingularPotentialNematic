%Finds if minima exist and where they occur for the MaierSaupe potential
%using the Ball-Majumdar Field theory
function [Iso,f0,Nem,Sn,fS] = MSMin(alph,LEB)
    GOLD = (sqrt(5)-1)/2;
    P = 0;
    %Find whether or not a min exists at 0
    fa = MaierSaupe(-0.05,P,alph,LEB);
    fb = MaierSaupe(0.05,P,alph,LEB);
    f0 = MaierSaupe(0,P,alph,LEB);
    if (f0 < fa) && (f0 < fb)
        Iso = true;
    else
        Iso = false;
    end
    %Find whether or not a min exists at a larger S value
    a = 0.2;
    c = 0.995;
    b = a + GOLD*(c-a);
    fa = MaierSaupe(a,P,alph,LEB);
    fc = MaierSaupe(c,P,alph,LEB);
    fb = MaierSaupe(b,P,alph,LEB);
    %try to bracket min in a couple steps
    if (fb < fa) && (fb < fc)
        Nem = true;
    elseif (fb > fa)
        c = b - (1-GOLD)*(b-a);
        fc = MaierSaupe(c,P,alph,LEB);
        [b,c,fb,fc] = deal(c,b,fc,fb);
        if (fb < fa) && (fb < fc)
            Nem = true;
        else
            a = b + GOLD*(c-b);
            fa = MaierSaupe(a,P,alph,LEB);
            [b,a,fb,fa] = deal(a,b,fa,fb);
            if (fb < fa) && (fb < fc)
                Nem = true;
            else
                Nem = false;
                Sn = b;
                fS = fb;
            end
        end
    else
        a = b + GOLD*(c-b);
        fa = MaierSaupe(a,P,alph,LEB);
        [b,a,fb,fa] = deal(a,b,fa,fb);
        if (fb < fa) && (fb < fc)
            Nem = true;
        else
            c = b - (1 - GOLD)*(b-a);
            fc = MaierSaupe(c,P,alph,LEB);
            [b,c,fb,fc] = deal(c,b,fc,fb);
            if (fb < fa) && (fb < fc)
                Nem = true;
            else
                Nem = false;
                Sn = b;
                fS = fb;
            end
        end
    end
    %Find min if Nem=true
    if Nem
        R = 1-GOLD;
        x0 = a;
        x3 = c;
        if (c-b > b-a)
            x1 = b;
            x2 = b + R*(c-b);
        else
            x2 = b;
            x1 = b - R*(b-a);
        end
        f1 = MaierSaupe(x1,P,alph,LEB);
        f2 = MaierSaupe(x2,P,alph,LEB);
        while (x3 - x0 > 1e-7*(abs(f1) + abs(f2)))
            if (f2 < f1)
                [x0,x1,x2] = deal(x1,x2,GOLD*x2+R*x3);
                [f1,f2] = deal(f2,MaierSaupe(x2,P,alph,LEB));
            else
                [x3,x2,x1] = deal(x2,x1,GOLD*x1+R*x0);
                [f2,f1] = deal(f1,MaierSaupe(x1,P,alph,LEB));
            end
        end
        if (f1 < f2)
            Sn = x1;
            fS = f1;
        else
            Sn = x2;
            fS = f2;
        end
    end
    
end