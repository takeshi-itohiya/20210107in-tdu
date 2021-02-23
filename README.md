//# 20210107in-tdu
//CDMAシステムの上り回線をモンテカルロ法でのシミュレーションする。シミュレーションの結果はｇｎｕｐｌｏｔで見る。

//使用言語　Ｃ言語

#include<stdio.h>//入出力
#include<math.h>
#include<stdlib.h>//乱数および強制終了

#define PI 3.14159265359//円周率

int main(void){

//関数のプロトタイプ宣言
	double G(double);//アンテナの利得を求める関数
	double Okumura_Hata(double,double,double,double);//奥村秦式
	double box_mullerd_method(void );//ボックスミュラー法による正規乱数
	double pythagorean_theorem(double,double);//三平方の定理により 直角三角形の斜辺でない辺から斜辺を求める。
	
//	使用する変数
//	一定の値をとる定数
	const int NB=19;//基地局の数
	const double hb=60;//基地局の高さ[m]
	const double hm=1.5;//端末の高さ[m]
	const double f=2000;//周波数は200MHZ,8000MHZ,2GHzとする。

	const double micro= 0.4;//有音率　SIRの計算で使用
	const double pg=3480/12.2;//システムパラメータ　拡散後の帯域幅/拡散前の帯域幅SIRの計算で使用
	const double beta=2; //送信電力増加率SIRの計算で使用
	const double PN_log= -129;//雑音電力[dBm]
	const double PN=pow(10,PN_log/10);//雑音電力を真数にしたもの
	const double SIR_req=6;//所要SIR[dB]
	const double Pmax=30;//最大送信電力[dBm] 送信電力制御で使用
	const double Pmin=-50;//最小送信電力[dBm]
	const int user_number_max=NB*200;//ユーザー数の最大値 基地局の数と１セル当たりのユーザ数最大値の積
	const double R=1.5;//セル半径[km]
	const double interstation_distance=2*R*cos(PI/6);//局間距離(リングの半径でもある。)[km]は、2*R*cos(π/6)である。π/6=30[deg]
	const double sigma=6;//短区間変動の標準偏差[dB] 6dB
	const double theta0=atan((hb-hm)/R)*180/PI;//セル端見込み角[deg]
	const double alpha=theta0+1;//チルト角(固定)[deg]
	const double epsilon=0.001;//SIRの許容誤差[dB]

	//計算の過程において変更してはいけないが、constでは扱えない変数
	double base_x[NB];//基地局座標x
	double base_y[NB];//基地局座標y


// 変化する変数
	//int N;//1セル当たりのユーザ 1から200まで
	int user_number;//ユーザ数
	int Nita;//イタレーション回数
	double Pt_log[user_number_max];//送信電力のログ[dBm] 　ユーザ毎の接続先の決定に使用 　初期値は0[dBm]
	double Pt[user_number_max];//送信電力
	double x[user_number_max];//ユーザのx座標[km]
	double y[user_number_max];//ユーザのY座標[km]
	int counter_centerbase;//中央の基地局に接続する端末の数
	int counter[NB];//それぞれの基地局に接続する端末の数
	int access_point[user_number_max];//接続先の場所 0が中央。
	double temp_Pr_log[user_number_max];//受信電力のログ[dBm]
	double Pr_log[user_number_max];//最大受信電力のログ
	double Pr[user_number_max];//最大受信電力の真数
	double temp_Gain_log[user_number_max];//ユーザのゲインのログ [dB]
	double Gain_log[user_number_max];//受信電力が最大時のユーザのゲインのログ
	double temp_Loss_log[user_number_max];//ユーザの伝送損失の対数[dB]
	double Loss_log[user_number_max];//受信電力が最大時の伝送損失のログ
	double temp_X_log[user_number_max];//短区間変動のログ
	double X_log[user_number_max];//受信電力最大時の短区間変動のログ
	double d;//基地局と ユーザの距離
	double theta;//角度
	double SIR[user_number_max];//各ユーザのSIR
	double SIR_log[user_number_max];//各ユーザのSIRの対数[dB]
	double delta_si[user_number_max];//各ユーザのSIRと所要
	int state_x;//所要値を満たすユーザの数
	int state_y;//品質劣化のユーザ数
	int state_z;//品質過剰のユーザ 数
	double temp_ratio;//品質劣化率の値
	double sum_ratio[200];//品質劣化率の合計
	double ratio[200];//品質劣化率(の平均)
	double sum_Pr[200];//すべてのユーザの受信電力の合計
	double temp_Pr_ave[200];//平均受信電力
	double sum_Pr_ave[200];//イタレーション回数に対する平均受信電力の
	double Pr_ave[200];//ユーザ数に対する平均受信電力
	double Pr_ave_log[200];//平均受信電力[dBm]
	
	
//ファイルポインタ
	FILE *fp_test1;//確認用のデータ
	FILE *fp_test2;//確認用のデータ
	FILE *fp_result;//結果のしゅつりょく　

//変数の実行
	
	{//ブロック１　基地局座標の決定

	//ブロック内の変数宣言
	
		int ring_number=1;//注目するリング//基地局座標の決定に使う。
		int a1=1;//注目するリングを移したときの最後の数 初期値は1とする。
		double theta_putbase=2*PI/3;//基地局を配置する際に使う角度 初期値120
		FILE *fp_test_base;//基地局のテストをするファイル

		base_x[0]=0;//中心基地局の座標は0。
		base_y[0]=0;
		//一つ目の基地局が
		base_x[1]=interstation_distance;//
		base_y[1]=0;//

		for(int i1=2;i1<NB;i1++){//基地局の番号
			if(i1-a1==ring_number*6){//リング内の番号-1が注目するリングの番号の６倍で割り切れたとき(六角形を一周したとき)
 				theta_putbase=2*PI/3;//角度を2/3π(120[deg])にする。
				ring_number+=1;//注目するリングの番号をインクリメント
				a1=i1;//注目する番号をリングが移動した最後の番号とする。
				base_x[i1]=(double)ring_number*interstation_distance;//x座標は局間距離の注目するリングの番号倍 とする。位置は、右端となる。
				base_y[i1] =0;//
				//printf("i1=%d\tx=%f\ty=%f\tringnumber=%d\ta1=%d\t一周\n",i1,base_x[i1],base_y[i1],ring_number,a1);//テスト用

			}else{
		 		base_x[i1]=(double)interstation_distance*cos(theta_putbase)+base_x[i1-1];//特定の角度に局間距離だけ動かす。
				base_y[i1]=(double)interstation_distance*sin(theta_putbase)+base_y[i1-1];//
				//printf("i1=%d\tx=%f\ty=%f\tringnumber=%d\ta1=%d\n",i1,base_x[i1],base_y[i1],ring_number,a1);
				//printf("theta=%f\n",theta_putbase*180/PI);
				if((i1-a1)%ring_number==0){//リング内の番号-1が注目するリングの番号で割り来れるとき(六角形の角)(最初は六角形の角にくるはず。
					theta_putbase+=(2*PI/6);//角度60[deg]足す
				}	
			}
			//
		}
		fp_test_base=fopen("test_base.txt","w");//課題１のファイルを開く
		if(fp_test_base==NULL){//ファイルが開けないとき
			printf("基地局位置のテストをオープンできませんでした。\n");
			return 1;
		}else{
		printf("基地局位置のテストファイルをオープンしました。\n");
		}

		for(int i1=0;i1<NB;i1++){
			fprintf(fp_test_base,"%f\t%f\n",base_x[i1],base_y[i1]);//
		}

		fclose(fp_test_base);//ファイルを閉じる。
		 printf("基地局位置のテストファイルをクローズしました。\n");
	}//ブロック１終了




	//ユーザを配置して計算する。
 	for(int N=1;N<=200;N++){//ユーザ数Nは1から200まで変化させる。200までは計算する。
		//printf("N=%dの計算開始",N);//実行の速度の確認
		sum_ratio[N-1]=0;//イタレーションに対する 品質劣化率の合計
		sum_Pr_ave[N-1]=0;//一回のイタレーションの平均受信電力の合計の初期化

		for(Nita=1;Nita<11;Nita++){//イタレーション回数は10
			do{//中心基地局に接続するユーザが0ならここまで戻る。
				counter_centerbase=0;//中心基地局に接続するユーザ数を初期化
				for(int i1;i1<user_number_max;i1++){//配列の初期化
					Pt_log[i1]=0;
					Pt[i1]=1;
				}
		
				for(int usernumber=0;usernumber<N*NB;usernumber++){//各ユーザの座標を決定
					do{
						x[usernumber]=(double)rand()/RAND_MAX*12-6;//±6との一様乱数とする
						y[usernumber]=(double)rand()/RAND_MAX*12-6;//±6との一様乱数とする
						//printf("i=%dztx=%f\ty=%f\n",usernumber,x[usernumber],y[usernumber]);//テスト用
					}while( (x[usernumber]*x[usernumber]+y[usernumber]*y[usernumber]) >36);//半径6km以内(xの２乗とyの２乗の和が36を超える)に入っていない場合はもう一度処理をやり直す。
					//printf("i=%dztx=%f\ty=%f\n",usernumber,x[usernumber],y[usernumber]);//テスト用
				}

				if(N==1&&Nita==1){//ユーザ数１のときユーザの位置を確認する
						FILE *fp_test_userN1;//ファイルポインタ
						fp_test_userN1=fopen("test_userN1.txt","w");//課題１のファイルを開く
					if(fp_test_userN1==NULL){//ファイルが開けないとき
						printf("N=1のときのユーザ配置のテストをオープンできませんでした。\n");
						return 1;
					}else{
						printf("N=1のときのユーザ配置のテストファイルをオープンしました。\n");
					}
					for(int usernumber=0;usernumber<N*NB;usernumber++) {
						fprintf(fp_test_userN1,"%f\t%f\n",x[usernumber],y[usernumber]);//ユーザごとの位置を記録する。
					}
					fclose(fp_test_userN1);
					printf("N=1のときのユーザ配置のテストファイルをクローズしました。\n");//
				}else if(N==15&&Nita==1){//ユーザ数が15
					FILE *fp_test_userN15;//ファイルポインタ
					fp_test_userN15=fopen("test_userN15.txt","w");//課題１のファイルを開く
					if(fp_test_userN15==NULL){//ファイルが開けないとき
						printf("N=15のときのユーザ配置のテストをオープンできませんでした。\n");
						return 1;
					}else{
						printf("N=1のときのユーザ配置のテストファイルをオープンしました。\n");
					}
					for(int usernumber=0;usernumber<N*NB;usernumber++) {
						fprintf(fp_test_userN15,"%f\t%f\n",x[usernumber],y[usernumber]);//ユーザごとの位置を記録する。
					}
					fclose(fp_test_userN15);
					printf("N=15のときのユーザ配置のテストファイルをクローズしました。\n");//
				}else if(N==200&&Nita==1){//ユーザ数が200のとき
					FILE *fp_test_userN200;//ファイルポインタ
					fp_test_userN200=fopen("test_userN200.txt","w");//課題１のファイルを開く
					if(fp_test_userN200==NULL){//ファイルが開けないとき
						printf("N=200のときのユーザ配置のテストをオープンできませんでした。\n");
						return 1;
					}else{
						printf("N=200のときのユーザ配置のテストファイルをオープンしました。\n");
					}
					for(int usernumber=0;usernumber<N*NB;usernumber++) {
						fprintf(fp_test_userN200,"%f\t%f\n",x[usernumber],y[usernumber]);//ユーザごとの位置を記録する。
					}
					fclose(fp_test_userN200);
					printf("N=200のときのユーザ配置のテストファイルをクローズしました。\n");//
				}

				for(int i1=0;i1<NB*N;i1++){//各ユーザ毎に
					for(int i2=0;i2<NB;i2++){//基地局毎に角ユーザの受信電力を計算し、最も近い基地局を求める。
						d=pythagorean_theorem(x[i1]-base_x[i2],y[i1]-base_y[i2]);//三平方の定理の関数を使用し、距離を計測
						theta=atan((hb-hm)/d)*180/PI;//距離と携帯電話の基地局の高さより 角度の計算[deg]
						temp_X_log[i1]=box_mullerd_method()*sigma;//短区間変動の計算
						temp_Gain_log[i1]=G(theta-alpha);//角度よりアンテナ利得の計算
						temp_Loss_log[i1]=Okumura_Hata(d,f,hb,hm);//伝送損失の計算[dB]
						temp_Pr_log[i1]=Pt_log[i1]+temp_Gain_log[i1]-temp_Loss_log[i1]+temp_X_log[i1];//受信電力の計算
						if(i2==0 ||temp_Pr_log[i1]>Pr_log[i1]){//最初もしくは受信電力が受信電力の最大値を上回っていた場合各種最大値の設定
							Pr_log[i1]=temp_Pr_log[i1];//電力　
							Gain_log[i1]=temp_Gain_log[i1];//利得
							Loss_log[i1]=temp_Loss_log[i1];//伝送損失
							X_log[i1]=temp_X_log[i1];//短区間変動
							Pr[i1]=pow(10,Pr_log[i1]/10);//受信電力の対数を真数に変変換
							access_point[i1]=i2;//接続先を最も電力の大きいところとする。
						}
					}
					//printf("i1=%dのときのアクセスポイント%d\n",i1,access_point[i1]);//テスト用
					if(access_point[i1]==0){//接続先が中央の基地局なら
						counter_centerbase++;//中央の基地局に接続するユーザを1増やす。
					}

					//printf("i1=%dのときの中央の基地局に接続するユーザ=%d\n",i1,counter_centerbase);//テスト用
					counter[ access_point[i1] ]++;//つながっている先の基地局の番号のカウンタをインクリメント
				}

				//printf("N=%d\tcounter=%d\n",N,counter_centerbase);//テスト用
			}while(counter_centerbase==0);//中央に接続するユーザが0なら前に戻る

			if(Nita==1&&(N==1||N==15||N==200)){
				FILE *fp_test_center;//中央につながっている数
				if(N==1){
					fp_test_center=fopen("test_centerN1.txt","w");

				}else if(N==15){
					fp_test_center=fopen("test_centerN15.txt","w");
				}else{
					fp_test_center=fopen("test_centerN200.txt","w");
				}

				if(fp_test_center==NULL){//ファイルが開けないとき
					printf("中央に接続するユーザのテストをオープンできませんでした。\n");
					return 1;
				}else{
				printf("中央に接続するユーザのテストをオープンしました。\n");
				}

				for(int i1=0;i1<N*NB;i1++){
					if(access_point[i1]==0) fprintf(fp_test_center,"%f\t%f\n",x[i1],y[i1]);//接続先が中央ならユーザの座標を書き込む
				}
				fclose(fp_test_center);//ファイルを閉じる
				printf("中央に接続するユーザのテストをクローズしました。\n");

			}
			do{//電力制御
				double sum_base[NB];//基地局ごとの（送信電力*利得*短区間変動/伝送損失）合計
				state_x=0; state_y=0; state_z=0; //初期化
				for(int i1=0;i1<N*NB;i1++){//基地局ごとの初期化
					sum_base[i1]=0;
				}
				for(int i2=0;i2<N*NB;i2++){//各ユーザ毎に基地局ごとの(受信電力（送信電力*利得*短区間変動/伝送損失）を求めるための計算
					int a2;//ユーザの接続先
					sum_base[ access_point[i2] ]+=Pr[i2];//注目するユーザに接続する基地局にユーザの受信電力を足す。 
				}
				//各ユーザのSIRの計算

				
				for(int i3=0;i3<N*NB;i3++){//i3はSIRを計算するユーザお番号
					double sum_other_b=0;///注目する基地局意外の他の基地局の（送信電力*利得*短区間変動/伝送損失）の合計
					for(int i4=0;i4<NB;i4++){//各基地局ごとの(送信電力*利得*短区間変動/伝送損失を計算する。)の合計
						if( i4!=access_point[i3]){//注目するユーザの基地局でないとき
							sum_other_b+=sum_base[i4];//
						}

						SIR[i3]=(Pr[i3])/( PN+micro/pg* ( (sum_base[ access_point[i3] ]-Pr[i3])+beta*sum_other_b));// SIRの計算分母は注目するユーザの受信電力、分子

						if(SIR[i3]<=0){//もしSIRが0かマイナスと出た場合  プログラムを終了する。
							printf("i3=%dのときSIRの数値が0かマイナスと出たので計算を終了します。\n",i3);
							exit(1);
						}
						SIR_log[i3]=10*log10(SIR[i3]);//
					}
					delta_si[i3]=SIR_log[i3]-SIR_req;//SIRと所要値の差
					if(fabs(delta_si[i3])<epsilon){//SIRの誤差の絶対値が許容誤差を下回るとき
						state_x++;//
					}else if(Pt_log[i3]==Pmax && delta_si[i3]<0){//品質劣化のとき  送信電力最大かつSIRが所要値を下回る
						state_y++;
					}else if(Pt_log[i3]==Pmin && delta_si[i3]>0){//品質過剰のとき　送信ん電力最小かつSIRが所要値を上回る。
						state_z++;
					}
				}
				if((N==2||N==4)&&Nita==1){//テスト用コード
					printf("N=%d\t%d\t%d\tx=%d\ty=%d\tz=%d\n",N,N*NB-state_x-state_y-state_z ,state_x+state_y+state_z ,state_x,state_y,state_z);
					for(int i7=0;i7<N*NB;i7++){
						printf("Pt[%d]=%f[dBm] SIR=%f[dB] Pr[%d]=%f[dB]\t",i7,Pt_log[i7],SIR_log[i7],i7,Pr_log[i7]);//テスト用
					}
					printf("\n");
				}
				if(N*NB-state_x-state_y-state_z==0){//収束条件 すべてのユーザが所要値を満たすか、品質過剰か品質劣化であること
					temp_ratio=(double)state_y/(double)(state_x+state_y+state_z);//品質劣化率
					sum_Pr[N-1]=0;//変数の初期化
					sum_ratio[N-1]+=temp_ratio;//Nはユーザ数よって、-1を引いて0オリジンに合わせる。
					for(int i5=0;i5<NB*N;i5++){ //すべてのユーザの受信電力の合計
						sum_Pr[N-1]+=Pr[i5];//
					}
					temp_Pr_ave[N-1]=sum_Pr[N-1] /(N*NB);//平均受信電力の計算
					sum_Pr_ave[N-1]+=temp_Pr_ave[N-1];//
					break;//ループを脱出する。
				}else{//収束しないとき電力を調整
					for(int i6=0;i6<=N*NB;i6++){//各ユーザの送信電力を調整する。受信電力の計算も同時に行う。
						Pt_log[i6]-=delta_si[i6];//送信電力の対数値からSIRの所要値との誤差を引く.
						if(Pt_log[i6]>Pmax){//送信電力が最大値を超えるなら、最大値にする。
							Pt_log[i6]=Pmax;//
						}else if(Pt_log[i6]<Pmin){//送信電力が最小値未満なら最小値にする。
							Pt_log [i6]=Pmin;//
						}
						Pr_log[i6]=Pt_log[i6]+Gain_log[i6]+X_log[i6]-Loss_log[i6];//受信電力の対数の計算
						Pr[i6]=pow(10,Pr_log[i6]/10);//対数を真数に直す
					}
				}
			}while(1);//無限ループ



		}//ここまできたら、10回繰り替える。(イタレーション)
		ratio[N-1]=sum_ratio[N-1]/10;//イタレーション回数に対する
		Pr_ave[N-1]=sum_Pr_ave[N-1]/10;//イタレーション回数に対する平均の計算。
	}

	{ //このブロックで結果を出力する。
		FILE *fp_pr;//受信電力を出力
		FILE *fp_ratio;//品質劣化率を出力
		double Pr_ave_min=Pr_ave_log[0];//平均受信電力最大値 軸の設定に使う//テスト用
		double Pr_ave_max=Pr_ave_log[0];// 最小値
		double ratio_min=ratio[0];//品質劣化率の最小値
		double ratio_max=ratio[0];//品質劣化率の最大値
		double work;
		fp_pr=fopen("kadai6_pr.txt","w");//平均受信電力の出力ファイル　課題6のグラフとなる
		fp_ratio=fopen("kadai6_ratio.txt","w");//品質劣化率の出力ファイル 

		for(int i7=0;i7<200;i7++){//i7はユーザ数Nに１を引いた数と同じ
			if(Pr_ave[i7]<=0){//平均受信電力が0以下なら
				printf("ユーザ数%dのときの受信電力が不正な値です。\n",i7+1);//
				exit(1);//
			}
			Pr_ave_log[i7]=10*log10(Pr_ave[i7]);//平均受信電力の対数を取る。
			fprintf(fp_pr,"%d\t%f\n",i7+1,Pr_ave_log[i7]);//ユーザ数と平均受信電力を出力
			fprintf(fp_ratio,"%d\t%f\n",i7+1,1-ratio[i7]);//1-品質劣化率を出力
			//以下は軸の範囲の確認用
			if(Pr_ave_min>Pr_ave_log[i7]){
				 Pr_ave_min=Pr_ave_log[i7];
			}else if(Pr_ave_max<Pr_ave_log[i7]){
				Pr_ave_max=Pr_ave_log[i7];
			}
			if(ratio_min>ratio[i7]){
				ratio_min=ratio[i7] ;
			}else if(ratio_max<ratio[i7]){
				ratio_max=ratio[i7];
			}//ここまで軸の範囲の確認用
		}
		
		fclose(fp_pr);
		fclose(fp_ratio);
		//以下2行は軸の範囲の確認用
		printf("Prの最大値=%f\t最小値=%f\n",Pr_ave_max,Pr_ave_min);
		printf("ratioの最大値=%f\t最小値=%f\n",1-ratio_max,1-ratio_min);
	}

	return 0;//
}
 //----------------------------------------------------
double pythagorean_theorem(double x,double y){//ピタゴラスの定理を使用し入力の二乗和の平方根を求める。
	double z=sqrt(x*x+y*y);//	
	return z;//
}

double G(double phi){ //角度[deg]に対してアンテナの指向性を[dB]で返す関数。
	const double derad=PI/180;//デグリによる角度をラジアンに変換する。

	const int M=16;//M=16
	double x;//出力する数
	double phi_rad;//ラジアンに変換済みの角度

	if(fabs(phi)>40){//角度が40°以上のとき、 
		phi_rad=40*derad;//角度を40°とする。 
		//printf("alpha>40deg\n");//テスト用
	}else{//角度が40°以下なら、角度は
		phi_rad=phi*derad;//角度はそのまま
	}
	
	if(fabs(phi_rad)<asin(1/ ((double )M*2)) ){//sin^-1(1/2M)より角度が小さいとき
		if(phi_rad==0) phi_rad=  1e-10;//角度0を10^-10とする。（零除算を防ぐため)
		x=fabs( sin( (double)M*PI*sin(phi_rad))/sin(PI*sin(phi_rad)) ) * cos(PI/2*sin(phi_rad)) /cos(phi_rad)  / (double)M ;
	}else{  //角度が大きいとき(角度が40°以上の場合、角度は40°と設定されるので、90°のときの考慮は不要
		x=fabs( 1 / sin(PI*sin(phi_rad)) ) *cos(PI/2*sin(phi_rad)) /cos(phi_rad)  / (double)M ;
	}
//	printf("真数xは%fです。\n",x);//テスト用
	x=20*log10(x);//ｘをログに直す[dB]
	//printf("G=%f\n",x);テスト用
	return x;//
}

double Okumura_Hata(double D,double f,double hb,double hm){//奥村秦式  引数は(水平距離D[km],周波数f[MHz],基地局アンテナ高さhb[m]、受信局高さhf[m])
	double a;//
	double loss;//損失[dB]

	//printf("距離Dは%fです。fは%fです。hbは%fです。hmは%fです。\n",D,f,hb,hm);//
	if(D==0){///距離0のときはログの計算ができないので他の式を用いる。
		loss=20*log10(hb-hm)+20*log10(f)+20*log10(4*PI/3)-40;   //自由空間の伝搬の式を用いる。   
//伝搬式は20log10(d)+20log10(f)+20log10(4π/c) ただしdの単位はm fの単位は[Hz] cは光速3*10^-8[m/s]
//20*log10(基地局高さ-受信局高さ)+20*log10(f)+120+20*log10(4*PI/3)-160 より 20*log10(D)+20*log10(f)+20*log10(4*PI/3)+20


	}else {
		a=(1.1*log10(f)-0.7)*hm-(1.56*log10(f)-0.8);//中小都市の値を使う。
		loss=69.55+26.16*log10(f)-13.82*log10(hb)+(44.9-6.55 *log10(hb))*log10(D)-a;
	}
	//printf("損失は%fです。\n",loss);  //テスト用
	return loss;//出力
		
}

double box_mullerd_method(void){//ボックスミュラー法での正規乱数
	double x=(double) rand()/ RAND_MAX ; //[0,1]までの一様乱数
	double y= (double)rand()/ RAND_MAX ;//[0,1]までの一様乱数
	double out;//出力
//	printf("x=%f\ty=%f\n",x,y);//テスト用
	if(x>1 || x<0 || y>1 || y<0){//0,1および、異常な値が出たとき
		exit(1);//プログラム強制終了
	}else if(x==0||x==1){ //xが0か１のとき
		out=0;//0を出力する。
	}else{//
		out=sqrt(-2*log(x))*cos(2*PI*y);//
	}
	return out;
}
