function ddy_num = Func_integration_function_drone(dy, Wrench, Mt_num, C_num)

% Aa_num = A_num(:, 1:2);
% Ap_num = A_num(:, 3:4);
% dAa_num = dA_num(:, 1:2);
% dAp_num = dA_num(:, 3:4);

% Mt_aa_num = Mt_num(1:2,1:2);
% Mt_ap_num = Mt_num(1:2,3:4);
% Mt_pa_num = Mt_num(3:4,1:2);
% Mt_pp_num = Mt_num(3:4,3:4);
% C_a_num = C_num(1:2);
% C_p_num = C_num(3:4);
% Mt_lambda_num = [Mt_pp_num, -Ap_num.'   ;
%                  Ap_num,     zeros(2,2)];
% dya_num = [dy(1); dy(2)];
% Active_Mt_num = Mt_aa_num - [Mt_ap_num, -Aa_num.']*(Mt_lambda_num\[Mt_pa_num; Aa_num]);
% Active_C_num = C_a_num - [Mt_ap_num, -Aa_num.']*(Mt_lambda_num\[C_p_num; (dAa_num-dAp_num*(Ap_num\Aa_num))*dya_num]);
% ddya_num = Active_Mt_num\(F_b-Active_C_num);
% ddyp_num = -Mt_lambda_num\[Mt_pa_num*ddya_num+C_a_num; Aa_num*ddya_num+(dAa_num-dAp_num*(Ap_num\Aa_num))*dya_num];
% ddyp_num = ddyp_num(1:2);
% ddy_num = [ddya_num; ddyp_num];

Active_Mt_num = Mt_num;
Active_C_num = C_num;
ddy_num = Active_Mt_num\(Wrench-Active_C_num);

end
