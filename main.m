clearvars; clc;

R  = 800;
L = 50;

A = Klotoid.findA(800, 50);

e = Klotoid(R, L, A);

[t, T, alpha, b, BS] = deltaTegetUzunlugu(e, 60);


%%
R = 200;
delta = 80;
[To, L] = Klotoid.delta2L(delta, R);
A = Klotoid.findA(R, L);

klotoid = Klotoid(R, L, A);

[T_, Z_, BS_] = kurbunTamamiKlotoid(klotoid);

%%

delta = 60;
R = 200;
A = 80;

L = Klotoid.findL(R, A);

e = Klotoid(R, L, A);

[t, T, alpha, b, BS] = deltaTegetUzunlugu(e, delta);

KB1 = 1081.15;
aralik = 20;
[KS1, B, KS2, KB2, app] = kilometraj(e, KB1, b, aralik);

anaParametreler = Klotoid.findMainParams(app, aralik);
