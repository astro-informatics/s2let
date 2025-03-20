% s2let_perf_tests

ell = [ 4 8 16 32 64 128 256 512 1024 2048 4096 ];
accuracy = [ ...
    2.5e-15 ... % 4
    3.97e-15 ... % 8
    8.0e-15 ... % 16
    2.1e-14 ... % 32
    4.10e-14 ... % 64
    7.32e-14 ... % 128
    1.38e-13 ... % 256
    2.31e-13 ... % 512
    5.2e-13 ... % 1024
    9.3e-13 ... % 2048
    2.0e-12 ... % 4096
    ];
speed = [ ...
    0.00019 ...  % 4
    0.00038 ... % 8
    0.00160 ... % 16
    0.0060 ... % 32
    0.047 ... % 64
    0.32 ... % 128
    2.5 ... % 256
    20.5 ... % 512
    195 ... % 1024
    1750 ... % 2048
    15000 ... % 4096
    ];
accuracy_multires = [ ...
    2.5e-15 ... % 4
    3.91e-15 ... % 8
    8.0e-15 ... % 16
    2.1e-14 ... % 32
    3.97e-14 ... % 64
    7.32e-14 ... % 128
    1.37e-13 ... % 256
    2.37e-13 ... % 512
    5.3e-13 ... % 1024
    9.3e-13 ... % 2048
    2.0e-12 ... % 4096
    ];
speed_multires = [ ...
    0.00021 ...  % 4
    0.00040 ... % 8
    0.00105 ... % 16
    0.00300 ... % 32
    0.018 ... % 64
    0.1 ... % 128
    0.75 ... % 256
    5.2 ... % 512
    45.0 ... % 1024
    400 ... % 2048
    3000 ... % 4096
    ];
nb_samples = ell;
nb_samples_pow = str2mat('t04', 't08', 't16', 't32', 't64', 't28', 't56', 't12', 't24', 't48', 't96');

figure('Position',[100 100 700 600])

subplot(2,1,1)
loglog( nb_samples, (1e-15)*ell, 'red','LineWidth', 2 )
hold on
loglog( nb_samples, accuracy, '--oblack', 'LineWidth', 2, 'MarkerSize',10)
loglog( nb_samples, accuracy_multires, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
xlabel('x1','FontSize',20)
ylabel('y1','FontSize',20)
set(gca,'XTick',nb_samples,'XTickLabel',nb_samples_pow,'FontSize',20)
set(gca,'YTick',[1e-16 1e-15 1e-14 1e-13 1e-12 1e-11] ,'FontSize',20)
axis([3e0 5300 1e-15 1e-11])
grid minor
%title('Accuracy of overall transform','FontSize',20)

subplot(2,1,2)
loglog( nb_samples, (1e-6)*ell.^3, 'red','LineWidth', 2 )
hold on
loglog( nb_samples, speed, '--oblack', 'LineWidth', 2, 'MarkerSize',10)
loglog( nb_samples, speed_multires, '-oblack', 'LineWidth', 2, 'MarkerSize',10)
grid on
xlabel('x2','FontSize',20)
ylabel('y2','FontSize',20)
axis([3e0 5300 1e-5 1e5])
set(gca,'YTick',[1e-6 1e-4 1e-2  1e0 1e2 1e4 1e6] ,'FontSize',20)
set(gca,'XTick',nb_samples,'XTickLabel',nb_samples_pow,'FontSize',20)
grid minor
%title('Speed of overall transform','FontSize',20)

