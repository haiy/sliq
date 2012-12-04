/*
-----------------------------------------------------
Author :haiyangfu 
Date:	Tue 14 August 2012
Description:
	A sliq classification program.
Revised:
	1.Gini Index calculate method
	2.histogram[][][]
	3.stack part Aug 8
	4.Gini_Table Aug 9
	5.Depth first and breadth first traverse
	6.father link and  attr index Aug 13
	7.matthews correlation and wighted accurary 13
	8.done 14 Aug
Compile: gcc -o sliq sliq_final.c -lm
Usage  : ./sliq  path/to/data/file
----------------------------------------------------
*/

#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<stdlib.h>
#include<math.h>

/*-----------------------------------------*/
#define log2(x)  log10(x)/log10(2)
/*******************************************/
//data file types definition
#define CacheSize 20000
typedef float DataType;
typedef struct AttrTable{
	DataType AttrValue;
	int  Seq;
}AttrTable;
AttrTable * Attr_Table;
FILE * DATA_FILE;
char * LineCache;
DataType *Compound;
DataType *Attr;
int AttrCount;
int CompCount;
/*******************************************/
//tree node data types
/*------Histogram--------*/
#define L_line 0
#define R_line 1
#define Active_row 0
#define Inactive_row 1
/*--------Label----------*/
#define Active 1
#define Inactive -1
typedef struct ClassTable{
	int NodePosition;
	DataType Label;
}ClassTable;
typedef ClassTable * ClassTablePtr;
typedef struct TreeNode{
	DataType BestSplit;
	float GiniIndex;
	float acc;    //Classification accuricy
	float mcc;    //Matthews Coefficient
	float mic;    //mutal information coefficient
	float wmcc;   //weighted MCC
	float wacc;   //weighted ACC
	float wmic;   //weighted MIC
	int FNodePos;
	int AttrIndex;
	int Histogram[2][2][2];	//histogram[0], the temporary and histogram[1]the best.
	int NodePosition;
	int PureNode;
	struct TreeNode *lchild,* rchild ,* father ;
}TreeNode;
typedef TreeNode * TreeNodePtr;
typedef TreeNodePtr ElemType;
int NodeCount=0;
int LayerCount=0;
ClassTablePtr  Class_Table;
/*********************************************/
//queue node data types  definition
typedef struct QueueNode{
	ElemType QNode_Value;
	struct QueueNode * next;
}QNode,* QNodePtr;
typedef struct Queue{
	QNodePtr front;
	QNodePtr rear;
}Queue;
/*------------------ -----------------------*/
//stack data structure
typedef TreeNodePtr  StackType ;
typedef struct StackNode{
	StackType Stack_Value;
	struct StackNode * front;
}StackNode,* StackNodePtr;
typedef struct Stack{
	StackNodePtr Base;
	StackNodePtr Top;
}Stack;
/*-------------------------------------------*/
Stack  *S;
int StackSize=0;
/*-------------------------------------------*/
//gini table structure
typedef struct GiniTable{
	float GiniIndex;
	int Seq;
}GiniTable;
/*-------------------------------------------*/
GiniTable * Gini_Table;
int his[2][2];
/*-------------------------------------------*/
//RI table structure
typedef struct RI{
	int AttrIndex;//Attribute index
	int HIS[2][2];
	float MCC;    //Matthews Coefficient
	float MIC;    //mutal information coefficient
	float ACC;    //overall accuray
	float wMCC;   //weighted MCC
	float wACC;   //weighted ACC
	float wMIC;   //weighted MIC
}RI, * RITable;
RITable RI_Table;
/*-------------------------------------------*/
//tree methods
typedef TreeNodePtr StackType;
void Deep_Traverse(TreeNodePtr root,void (* Visit_Node)(TreeNodePtr) );
void View_Node(TreeNodePtr n);
void BFS_Traverse(Queue *Q ,TreeNodePtr root,void(* Visit_Node)(TreeNodePtr));
void (* Visit_Node)(TreeNodePtr);
/*-------------------------------------------*/

/* Functions Announcement */
/*-------------------------------------------*/
// Stack Methods
void Initial_Stack();
void Push_Stack(StackType Stack_Value);
StackType Out_Stack();
StackType Get_Front_Stack();
int Empty_Stack();
StackNodePtr Search_Stack(StackType Stack_Value);
void Traverse_Stack();
void Destroy_Stack();
/*------------------------------------------*/
//Calculatation Part 
void Get_Class_Table();
void Cal_Gini_Index(Queue *Q);
QNodePtr Search_Queue(Queue *Q,int NodeIndex);
void Sync_Histogram(int a[2][2],int b[2][2]);
float Gini_Index(int Histogram[2][2]);
void Update_Labels(Queue *Q);
void Purify_Queue(Queue * Q);
void Generate_Child_Node(Queue *Q);
void Initialize_Child_Node(TreeNodePtr Node);
void Add_Attr_Index(Queue *Q);
void MCC_ACC_Coefi();
/*-------------------------------------------*/
// Queue methods
void InitialQueue(Queue * Q);
void EnQueue( Queue * Q,ElemType QNode_Value);
int CalNode( Queue * Q);
ElemType Get_Front_Queue(Queue *Q);
int EmptyQueue(Queue *Q);
void DestroyQueue(Queue *Q);
void TraverseQueue(Queue *Q);
ElemType Check_Queue(Queue *Q,int QPosition);
/*------------------------------------------*/
//Data Prepare
void Initial_Data();
int Count_Attr();
int Count_Comp();
void Get_One_Compound(int CompIndex);
void Get_One_Attr(int AttrIndex);
void View_Value(int choice);
void Destroy_Data();
//-> Attr Table
void Get_Attr_Table(int AttrIndex);
void View_Attr_Table();
void Sort_Attr_Table(int count);
void Change_Position_Attr(int a,int b);
void Quick_Sort_Attr(int count);
void Sort_Attr(int low, int high);
int Partion_Attr(int low, int high);
//-> Class Table
void View_Class_Table();
void View_Histogram(int a[2][2]);
//-> Gini Table
void Get_Gini_Table();
void Initial_His();
void His_State();
float One_Attr_Gini();
float Attr_Gini(int his[2][2]);
void Sort_Gini_Table();
void View_Gini_Table();
/*-------------------------------------------*/
//Get RI Table
void Initial_RI_Table();
void Destroy_RI_Table();
void Initial_HIS(int h[2][2]);
void Get_RI_Table(Queue *Q ,TreeNodePtr root);
void Cal_Node_AM(TreeNodePtr node);
void Cal_RI_wAM(TreeNodePtr node);
void Cal_RI_tAM();
int Total_HIS(int h[2][2]);
void Add_His(int a[2][2],int b[2][2]);
float MCC_Cal(int a[2][2]);
float ACC_Cal(int a[2][2]);
void View_RI_Table();
void Destroy_RI_Table();
/*---------------The Main function------------*/
/*For MIC Use Sep 12 ,2012
void main(int arg ,char *argv[]){
*/
void Sliq_Tree(int arg,char *argv[]){
/* Part I: Initialize the program */

/* Initial Step 1: variables announcement and malloc memories  */
	Queue Q;
	TreeNodePtr root;
	int AttrIndex=1,i;
	DATA_FILE=fopen(argv[1],"r");
	if(!DATA_FILE){
		puts("File Open Error!");
		exit(0);	
	}
	Initial_Data();
	InitialQueue(&Q);
	Initial_Stack();
	Get_Gini_Table();//Get gini table
	View_Gini_Table();
/* Initial Step 2: root node initialization.*/
	root=(TreeNodePtr)malloc(sizeof(TreeNode));
	Initialize_Child_Node(root);
	Get_Attr_Table(Gini_Table[0].Seq); //the root attr,the second row of the data file.
	printf("Now the root attr!row : %d\n",Gini_Table[0].Seq);
	View_Attr_Table();
	NodeCount=root->NodePosition=1;    //root node position and the node count!Important!
	root->AttrIndex=0;
	root->FNodePos=0;
	Get_Class_Table();
	for(i=0;i<CompCount;i++){
		if(Class_Table[i].Label==Active)
			root->Histogram[0][1][0]+=1;
		else
			root->Histogram[0][1][1]+=1;
	}
	View_Histogram(root->Histogram[0]);
/* Part II :Calculate  */
	LayerCount=1;
	i=1;
	EnQueue(&Q,root);
	TraverseQueue(&Q);
	while(!EmptyQueue(&Q)){          //Only if the queue is empty then stop .
		puts("\nNew Layer Begain...");
		puts("-------------Cal Gini index --------------------");
		Cal_Gini_Index(&Q);      //Calaulate all the none pure nodes' gini index in the queue.
		puts("-------------Child Node bearing-----------------");
		Generate_Child_Node(&Q); //Generate and put all child nodes in queue.
		TraverseQueue(&Q);
		puts("-------------Update Labels----------------------");
		Update_Labels(&Q);       //Update the class table and initialize the children's histogram.
		View_Class_Table();
		TraverseQueue(&Q);
		Add_Attr_Index(&Q);
		puts("-------------Purify Queue-----------------------");
		Purify_Queue(&Q);        //Make sure that only the none pure nodes in the queue.  
		TraverseQueue(&Q);
		if(i<AttrCount){  	 //Get next attribute to calculate.
			Get_Attr_Table(Gini_Table[i].Seq);
			printf("Get the %d attr !\n",Gini_Table[i].Seq);
		}
		else
			break;
		LayerCount++;
		i++;
		printf("This is the %d time !",i);
	}
	LayerCount-=1;
	puts(">>>>>>>>Begin traverse ....");
	Initial_RI_Table();
	Get_RI_Table(&Q,root);
	BFS_Traverse(&Q,root,View_Node);
	DestroyQueue(&Q);
	Destroy_Data();
	Destroy_Stack();
	Destroy_RI_Table();
}
/* main function over. */
/*------------------------------------------------------------*/
/* Functions below are support functions */

/* I: Data Prepare */
void Initial_Data(){
	int i=0;
	puts("Begin Initializing...");
	LineCache=(char *)malloc((sizeof(char))*(CacheSize));
	AttrCount=Count_Attr();
	CompCount=Count_Comp();
	Compound=(DataType *)malloc(sizeof(DataType)*(AttrCount));
	Attr=(DataType *)malloc(sizeof(DataType)*(CompCount));
	Attr_Table=(AttrTable *)malloc(sizeof(AttrTable)*CompCount);
	Gini_Table=(GiniTable *)malloc(sizeof(GiniTable)*AttrCount);
	printf("%d comps and %d attrs.\n",CompCount,AttrCount);
}
int Count_Attr(){
	int AttrCount=0;
	char * temp;
	fseek(DATA_FILE,0,SEEK_SET);
	if(fgets(LineCache,CacheSize,DATA_FILE)){
		temp=LineCache;
		while(*temp!='\0'&&*temp!='\n'){
			strtod(temp,&temp);
			++AttrCount;
			temp++;	
		}
	}
	return AttrCount;
}
int Count_Comp(){
	int CompCount=0;
	fseek(DATA_FILE,0,SEEK_SET);
	while(fgets(LineCache,CacheSize,DATA_FILE))
		++CompCount;
	return CompCount;
}
void Get_One_Compound(int CompIndex){
	int i=0;
	char *temp;
	fseek(DATA_FILE,0,SEEK_SET);
	while(fgets(LineCache,CacheSize,DATA_FILE)){
//		printf("Now i %d,CompIndex %d\n",i,CompIndex);
		if(i++==CompIndex){
			i=0;
			temp=LineCache;
			while(*temp!='\n'&&*temp!='\0'){
				Compound[i]=strtod(temp,&temp);
//				printf("Get %d: %f.\n",i,Compound[i]);
				temp++;
				i++;
			}
//			printf("Get Compound %d success!\n",CompIndex);
			break;
		}
	}
}
void Get_One_Attr(int AttrIndex){
	int i;
	for(i=0;i<CompCount;i++){
		Get_One_Compound(i);
//		printf("Compound %f\n",Compound[AttrIndex]);
		Attr[i]=Compound[AttrIndex];
//		printf("Attr[%d]:%f,Compound[%d]:%f\n",i,Attr[i],AttrIndex,Compound[AttrIndex]);
	}
//	printf("Get Attr %d success!\n",AttrIndex);
}
void View_Value(int choice){   // choice = 1 to view the front 10 value of a Compound
//	printf("Choice is %d\n",choice);
	if(choice==1){
		for(choice=0;choice<10;choice++)
			printf("%f\t",Compound[choice]);
		printf("\n");
	}
	else if(choice==2){
		for(choice=0;choice<10;choice++)
			printf("%f\t",Attr[choice]);
		printf("\n");
	}
	else
		printf("Wrong Input!ONLY 1 or 2 will be accepted!\n");

}
void Get_Attr_Table(int AttrIndex){  //Get sorted attribute table.
	int i,count;
	Get_One_Attr(AttrIndex);
	for(i=0;i<CompCount;i++){
		Attr_Table[i].AttrValue=Attr[i];
//		printf("AttrIndex %d Get the %d ,Compound %f ,attr value %f\n",AttrIndex,i,Attr[i],Attr_Table[i].AttrValue);
		Attr_Table[i].Seq=i;
	}
	count=CompCount;
	Sort_Attr_Table(count);	
}
/* Quick sort Part */
void Sort_Attr_Table(int count){
	Quick_Sort_Attr(count);
}
void Change_Position_Attr(int a ,int b){
	DataType temp_value;
	int temp_seq;
	temp_value=Attr_Table[a].AttrValue;
	temp_seq=Attr_Table[a].Seq;
	Attr_Table[a].AttrValue=Attr_Table[b].AttrValue;
	Attr_Table[a].Seq=Attr_Table[b].Seq;
	Attr_Table[b].AttrValue=temp_value;
	Attr_Table[b].Seq=temp_seq;
}
void Quick_Sort_Attr(int count){
	Sort_Attr(0,count-1);
}
void Sort_Attr(int low,int high){
	int Middle;
	if(low<high){
		Middle=Partion_Attr(low,high);
		Sort_Attr(low,Middle-1);
		Sort_Attr(Middle+1,high);
	}
}
int Partion_Attr(int low,int high){
	DataType key=Attr_Table[low].AttrValue;
	while(low<high){
		while(low<high&&Attr_Table[high].AttrValue>=key) --high;
		Change_Position_Attr(low,high);
		while(low<high&&Attr_Table[low].AttrValue<=key) ++low;
		Change_Position_Attr(low,high);
		Attr_Table[low].AttrValue=key;
	}
	return low;
}
/* Quick sort part over. */
void View_Attr_Table(){
	int i=0;
	for(i=0;i<CompCount;i++){
		printf("Position:%d\tValue:%f label:%f \n",Attr_Table[i].Seq,Attr_Table[i].AttrValue,Class_Table[Attr_Table[i].Seq].Label);
	}
}
void Destroy_Data(){
	puts("Begin to clean memory!");
	free(Attr);
	free(LineCache);
	free(Compound);
	free(Attr_Table);
	free(Gini_Table);
	fclose(DATA_FILE);
	puts("Success!");
}
/*-------------------------------------------------------*/

/* II : Methods of the Queue */

void InitialQueue(Queue * Q){
	Q->front=Q->rear=(QNodePtr)malloc(sizeof(QNode));
	if(Q->front==NULL)
		printf("Queue Initial Error!\n");
	Q->front->next=NULL;
}
void EnQueue( Queue * Q,ElemType QNode_Value){
	QNodePtr TempNode;
	TempNode=(QNodePtr)malloc(sizeof(QNode));
	TempNode->QNode_Value=QNode_Value;
	TempNode->next=NULL;
	Q->rear->next=TempNode;
	Q->rear=TempNode;
}
int CalNode( Queue * Q){
	QNodePtr p;
	int Num=0;
	p=Q->front->next;
	while(p){
		++Num;
		p=p->next;
	}
	return Num;
}
ElemType Check_Queue(Queue *Q,int QPosition){// the first one from "0"
	QNodePtr p=Q->front->next;
	int Num=0,t=CalNode(Q);
	if(QPosition>t){
		printf("The maximum length of the Queue is %d!\n",t);
		return NULL;
	}
	else{
		while(Num<t){
			if(Num==QPosition)
				return p->QNode_Value;
			else{
				p=p->next;
				++Num;
			}
		}
	}
}
ElemType Get_Front_Queue( Queue *Q){
	ElemType temp_value=Q->front->next->QNode_Value;
	QNodePtr temp_node=Q->front;
	Q->front=temp_node->next;
	free(temp_node);
	return temp_value;
}
int EmptyQueue(Queue *Q ){
	if(Q->front==Q->rear)
		return 1;
	else 
		return 0;
}
void DestroyQueue(Queue *Q){
	QNodePtr temp_node=Q->front;
	printf("Destrory Queue!\n");
	while(!EmptyQueue(Q)){
		Q->front=temp_node->next;
		free(temp_node);
		temp_node=Q->front;
	}
}
/*-----------------------------------------------------------*/

/* III : Functions about the Tree Node*/

void Get_Class_Table(){
//Class table .
	int i=0;
	Get_One_Attr(0);
	if(Attr[i]!=1&&Attr[i]!=-1){
		puts("Data File Format Error!\n->The first row must be labels like -1 or 1.");
		printf("the value is :%f\n",Attr[i]);
		exit(0);
	}
	Class_Table=(ClassTablePtr)malloc(CompCount*sizeof(ClassTable));
	for(i=0;i<CompCount;i++){
		Class_Table[i].NodePosition=1;
		Class_Table[i].Label=Attr[i];
	}
}
void View_Class_Table(){
	int i=0;
//	printf("*****************Now Class Table*******************\n");
	for(i=0;i<CompCount;i++){
		printf("NodePosition: %d\tLabel:%f\n",Class_Table[i].NodePosition,Class_Table[i].Label);
	}
}
void Cal_Gini_Index(Queue *Q){//Need the Attr_Table first.
//traverse the attribute table and update all the nodes' Gini index in the queue.
	int i,j,L_new,temp_Histogram[2][2]={0},flag=1;
	float NewGiniIndex;
	QNodePtr temp_Qnode;
	TreeNodePtr temp_Tnode;
	for(i=0;i<CompCount-1;i++){
//		printf("Traverse the attr table.attr i:%d\n",i);
		temp_Qnode=Search_Queue(Q,Class_Table[Attr_Table[i].Seq].NodePosition);
		if(temp_Qnode==NULL){
			printf("Pure Node ,Pass!\n");
		}
		else{
//			printf("find the %d item of attr table's node in the queue!\n",Attr_Table[i].Seq);
			temp_Tnode=temp_Qnode->QNode_Value;
			Sync_Histogram(temp_Histogram,temp_Tnode->Histogram[0]);
			if(Class_Table[Attr_Table[i].Seq].Label==Active){
				temp_Histogram[0][0]+=1;
				temp_Histogram[1][0]-=1;
			}
			else{
				temp_Histogram[0][1]+=1;
				temp_Histogram[1][1]-=1;
			}
			flag=temp_Histogram[1][0]+temp_Histogram[1][1];
			Sync_Histogram(temp_Tnode->Histogram[0],temp_Histogram);
			if(flag&&(Gini_Index(temp_Histogram) < temp_Tnode->GiniIndex)){
				Sync_Histogram(temp_Tnode->Histogram[1],temp_Histogram);
				temp_Tnode->GiniIndex=Gini_Index(temp_Histogram);
				temp_Tnode->BestSplit=(Attr_Table[i].AttrValue+Attr_Table[i+1].AttrValue)/2.00;
			}
		}
		printf("Cal the i:%d ,Value[i]:%f,Value[i+1]:%f,Gini:%f \n\n",i,Attr_Table[i].AttrValue,Attr_Table[i+1].AttrValue,temp_Tnode->GiniIndex);
	}
}
QNodePtr Search_Queue(Queue *Q,int NodeIndex){
//Find the correspond Node in the Queue.
	QNodePtr temp_Qnode=Q->front->next;
	int Count=CalNode(Q);
	int i=1;
	while(i<=Count){
		if(temp_Qnode->QNode_Value->NodePosition==NodeIndex)
			return temp_Qnode;
		else
			temp_Qnode=temp_Qnode->next;
			i++;
	}
	return NULL; // not found in the queue return null.
}
void Sync_Histogram(int a[2][2],int b[2][2]){
//make array a the same as b.
	int i=0,j=0;
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			a[i][j]=b[i][j];
}
void View_Histogram(int a[2][2]){
	int *p=&a[0][0],i=0;
	printf("Histogram:\n");
	while(i<4){
		printf("Value:%d\t",*(p+i++));
		if(i==2)
			printf("\n");
	}
	printf("\n");
}
float Gini_Index(int Histogram[2][2]){
//calculate the Gini Index according to the value in Histogram
	float Gini_L,Gini_R,pA,pI,pL,pR,Gini_Index;
	int Total_1=Histogram[0][0]+Histogram[0][1],Total_2=Histogram[1][0]+Histogram[1][1];
//	View_Histogram(Histogram);
	pA=(float)(Histogram[0][0])/(float)(Total_1);
	pI=(float)(Histogram[0][1])/(float)(Total_1);
	Gini_L=1-(pA)*(pA)-(pI)*(pI);
//	printf("Left :%f\t%f\t%f\t",pA,pI,Gini_L);
	pA=(float)(Histogram[1][0])/(float)(Total_2);
	pI=(float)(Histogram[1][1])/(float)(Total_2);
	Gini_R=1-(pA)*(pA)-(pI)*(pI);
//	printf("Right :%f\t%f\t%f\t\n",pA,pI,Gini_R);
	pL=(float)Total_1/(float)(Total_1+Total_2);
	pR=(float)Total_2/(float)(Total_1+Total_2);
	Gini_Index=pL*Gini_L+pR*Gini_R;
//	printf("Gini_Index :%f\t%f\t--->%f\t\n",pL,pR,Gini_Index);
	return Gini_Index;
}
/* Update Class Labels of the class table */
void Update_Labels(Queue *Q){
	int i=0;
	QNodePtr temp_node;
	TreeNodePtr temp_Tnode;
	for(i=0;i<CompCount;i++){//After the child node generated,update labels and initial histogram.
		if((temp_node=Search_Queue(Q,Class_Table[Attr_Table[i].Seq].NodePosition))){
			temp_Tnode=temp_node->QNode_Value;
			if(Attr_Table[i].AttrValue < temp_node->QNode_Value->BestSplit){ //less than the split point to left part.
				Class_Table[Attr_Table[i].Seq].NodePosition=temp_node->QNode_Value->lchild->NodePosition;
				if(Class_Table[Attr_Table[i].Seq].Label==Active)
					temp_node->QNode_Value->lchild->Histogram[0][1][0]+=1;
				else
					temp_node->QNode_Value->lchild->Histogram[0][1][1]+=1;
				temp_node->QNode_Value->lchild->PureNode=(temp_Tnode->lchild->Histogram[0][1][0] < temp_Tnode->lchild->Histogram[0][1][1])?temp_Tnode->lchild->Histogram[0][1][0]:temp_Tnode->lchild->Histogram[0][1][1];
			}
			else{//right part of the split node
				Class_Table[Attr_Table[i].Seq].NodePosition=temp_node->QNode_Value->rchild->NodePosition;
				if(Class_Table[Attr_Table[i].Seq].Label==Active)
					temp_node->QNode_Value->rchild->Histogram[0][1][0]+=1;
				else
					temp_node->QNode_Value->rchild->Histogram[0][1][1]+=1;

				temp_node->QNode_Value->rchild->PureNode=temp_Tnode->rchild->Histogram[0][1][0] < temp_Tnode->rchild->Histogram[0][1][1]?temp_Tnode->rchild->Histogram[0][1][0]:temp_Tnode->rchild->Histogram[0][1][1];
			}
		}
		else
//			printf("Update labels failure . Not found the AttrValue:%f Attr Seq:%d in the Queue.\n ",Attr_Table[i].AttrValue,Attr_Table[i].Seq);
			;
	}
}
/* Split Node */
void Generate_Child_Node(Queue *Q){//this one should be  before update labels
	QNodePtr temp_Qnode=Q->front->next;
	int Count=CalNode(Q);
	int i=1;
	while(i<=Count){//Generate the children nodes and push them into the queue.
		temp_Qnode->QNode_Value->lchild=(TreeNodePtr)malloc(sizeof(TreeNode));
		NodeCount+=1;
		temp_Qnode->QNode_Value->lchild->NodePosition=NodeCount;
		Initialize_Child_Node(temp_Qnode->QNode_Value->lchild);
		EnQueue(Q,temp_Qnode->QNode_Value->lchild);
		temp_Qnode->QNode_Value->rchild=(TreeNodePtr)malloc(sizeof(TreeNode));
		NodeCount+=1;
		temp_Qnode->QNode_Value->rchild->NodePosition=NodeCount;
		Initialize_Child_Node(temp_Qnode->QNode_Value->rchild);
		EnQueue(Q,temp_Qnode->QNode_Value->rchild);
		temp_Qnode->QNode_Value->lchild->FNodePos=temp_Qnode->QNode_Value->NodePosition;
		temp_Qnode->QNode_Value->rchild->FNodePos=temp_Qnode->QNode_Value->NodePosition;
		temp_Qnode->QNode_Value->lchild->father=temp_Qnode->QNode_Value->rchild->father=temp_Qnode->QNode_Value;
		temp_Qnode=temp_Qnode->next;
		i++;
	}
}
void Initialize_Child_Node(TreeNodePtr Node){
	int i=0,*p=&(Node->Histogram[0][0][0]);
	Node->BestSplit=0.00;
	Node->PureNode=1;  //PureNode = 1 means a none pure node.
	Node->GiniIndex=1.00;
	Node->lchild=Node->rchild=NULL;
	Node->AttrIndex=0;
	Node->FNodePos=0;
	Node->acc=Node->mcc=Node->wmcc=Node->wacc=0.00;
	while(i<6){
		*p++=0;
		i++;
	}
}
/*-------------------Add AttrIndex to node---------------------*/
void Add_Attr_Index(Queue *Q){
	int Count=CalNode(Q);
	QNodePtr temp=Q->front->next;
	int p=temp->QNode_Value->AttrIndex;
	while(temp){
		if(!temp->QNode_Value->BestSplit){
			temp->QNode_Value->AttrIndex=p+1;
		}
		temp=temp->next;
	}

}
/*-------------------Add AttrIndex Over------------------------*/
/*-------------------------------------------------------------*/
void Purify_Queue(Queue *Q){
	QNodePtr temp1,temp2; 
	int Count=CalNode(Q);
	int i=1;
	temp1=Q->front;
	temp2=temp1->next;
	printf("%d nodes in the queue!\n",CalNode(Q));
	while(temp2){
		if(temp2->QNode_Value->BestSplit||(!temp2->QNode_Value->PureNode)||(!(temp2->QNode_Value->Histogram[0][1][0]&&temp2->QNode_Value->Histogram[0][1][0]))){
			if(temp2->next){
				temp1->next=temp2->next;
			}
			else{
				Q->rear=temp1;
				temp1->next=NULL;
			}
			printf("NodePositon %d out the Queue!\n",temp2->QNode_Value->NodePosition);
			free(temp2);
			temp2=temp1->next;
		}
		else{
			temp1=temp2;
			temp2=temp1->next;
			if(temp1&&temp1->next){
				temp1->next;
			}
			else if(temp1==NULL){
				Q->rear=temp2;
				break;
			}
		}
	}
}
void TraverseQueue(Queue *Q){
	QNodePtr temp_node = Q->front->next;
	int Count=CalNode(Q);
	int i=1;
	while(i<=Count){
		printf("\n--->\nTraversing the %d node of %d,\n BestSplit:%f,PureNode:%d ,\nNodePosition :%d,NodeCount %d,GiniIndex: %f\n",i,Count,temp_node->QNode_Value->BestSplit,temp_node->QNode_Value->PureNode,temp_node->QNode_Value->NodePosition,NodeCount,temp_node->QNode_Value->GiniIndex);
		puts("Temporary Histogram:");
		View_Histogram(temp_node->QNode_Value->Histogram[0]);
		puts("Best Split Point:");
		View_Histogram(temp_node->QNode_Value->Histogram[1]);
		temp_node=temp_node->next;
		i++;
	}
}
/*---------------------------------------------------------------*/

/* IV Part : Stack Functions */

void Initial_Stack(){
	S=(Stack *)malloc(sizeof(Stack));
	S->Top=(StackNodePtr)malloc(sizeof(StackNode));
	S->Base=S->Top->front=(StackNodePtr)malloc(sizeof(StackNode));
}
void Push_Stack(StackType Stack_Value){
	StackNodePtr temp=(StackNodePtr)malloc(sizeof(StackNode));
	StackSize+=1;
	temp->Stack_Value=Stack_Value;
	temp->front=S->Top->front;
	S->Top->front=temp;
}
StackType Out_Stack(){
	StackType p=S->Top->front->Stack_Value;
	StackNodePtr temp_node=S->Top->front;
	S->Top->front=S->Top->front->front;
	free(temp_node);
	return p;
}
StackType Get_Front_Stack(){
	return S->Top->front->Stack_Value;
}
StackNodePtr Search_Stack(StackType Stack_Value){
	StackNodePtr temp=S->Top->front;
	int i=0;
	while(temp!=S->Base){
		if(temp->Stack_Value=Stack_Value)
			return temp;
		else
			return NULL;
	}
}
/*
void Traverse_Stack(){
	StackNodePtr temp=S->Top->front;
	int i=0;
	while(temp!=S->Base){
		printf("Index: %d Value :%d\n",i++,temp->Stack_Value);
		temp=temp->front;
	}
}
*/
int Empty_Stack(){
	if(S->Top->front==S->Base)
		return 1;
	else
		return 0;
}
void Destroy_Stack(){
	while(!Empty_Stack()){
		Out_Stack();
	}
}
/*--------------------------------------------------------------*/

/* V Part : Gini Table */

void Get_Gini_Table(){ //the first attr is label.
	Get_Class_Table();
	int i;
	Initial_His();
	View_Histogram(his);
	for(i=0;i<AttrCount-1;i++){ // Attention, the first row will be ignored !
		Get_Attr_Table(i+1);
		Gini_Table[i].Seq=i;
		Gini_Table[i].GiniIndex=One_Attr_Gini();
	}
	View_Gini_Table();
	puts("sort gini table");
	Sort_Gini_Table();
}
void View_Gini_Table(){
	int i;
	for(i=0;i<AttrCount-1;i++){
		printf("Attr:%d Sequence:%d,Gini Index:%f\n",i,Gini_Table[i].Seq,Gini_Table[i].GiniIndex);
	}
}
void Sort_Gini_Table(){
	int i,count;
	Attr_Table=(AttrTable *)realloc(Attr_Table,sizeof(AttrTable)*AttrCount);
	if(Attr_Table==NULL)
		printf("Realloc failed!\n");
	for(i=0;i<AttrCount-1;i++){      // "AttrCount - 1 " !!!
		Attr_Table[i].Seq=Gini_Table[i].Seq;
		Attr_Table[i].AttrValue=Gini_Table[i].GiniIndex;
	}
	count=AttrCount-1;               // Here !
	Sort_Attr_Table(count);
	for(i=0;i<AttrCount-1;i++){      // "AttrCount - 1 " !!!
		Gini_Table[i].Seq=Attr_Table[i].Seq;
		Gini_Table[i].GiniIndex=Attr_Table[i].AttrValue;
	}
	Attr_Table=(AttrTable *)realloc(Attr_Table,sizeof(AttrTable)*CompCount);
}
void Initial_His(){
	int i;
	for(i=0;i<CompCount;i++){
		if(Class_Table[i].Label==Active)
			his[1][0]+=1;
		else
			his[1][1]+=1;
	}
}
void His_State(){
	his[1][0]=his[0][0]+his[1][0];
	his[1][1]=his[0][1]+his[1][1];
	his[0][0]=his[0][1]=0;
}
float One_Attr_Gini(){
	int i;float Gini=1;
	His_State();
	for(i=0;i<CompCount-1;i++){
		if(Class_Table[Attr_Table[i].Seq].Label==Active){
			his[0][0]+=1;
			his[1][0]-=1;
		}
		else{
			his[0][1]+=1;
			his[1][1]-=1;
		}
		if(Gini_Index(his) < Gini){//if Gini denotes Info Gain then here should be ">" 
			Gini=Gini_Index(his);
		}
	}
//	printf("Now the biggest one %f\n",Gini);
	return Gini;
}
/* it is too puzzle to use information gain ...
float Attr_Gini(int his[2][2]){
	float x=his[0][0],y=his[0][1],p=his[1][0],q=his[1][1];
	float T=x+y+p+q,T1=x+y,T2=T-T1;
	printf("%d,%d,%d,%d,%f\n",his[0][0],his[0][1],his[1][0],his[1][1],T);
	float Info_XT,Gain,Info_T,Info_1,Info_2;
	if(!(x&&y&&y&&p&&q))
		return 0;
	Info_T=-(((x+p)/T)/log2((x+p)/T)+((y+q)/T)/log2((y+q)/T));
	Info_1=-((x/T1)/log2(x/T1)+(y/T1)/log2(y/T1));
	Info_2=-((p/T2)/log2(p/T2)+(q/T2)/log2(q/T2));
	Info_XT=(T1/T)*Info_1+(T2/T)*Info_2;
	printf("Info_T:%f,Info_1:%f,Info_2:%f,Info_XT:%f\n",Info_T,Info_1,Info_2,Info_XT);
	printf("\nGain:%f\n",Info_T-Info_XT);
	return Gain=Info_T-Info_XT;
}
*/
/*--------------------Tree Traverse-----------------------------------*/
void Deep_Traverse(TreeNodePtr root,void (* Visit_Node) (TreeNodePtr)){   // left depth first .
	TreeNodePtr p=root;
	Push_Stack(root);
	printf("Now stack size %d !\n",StackSize);
	puts("root in..");
	while(!Empty_Stack()){
		while(p&&Get_Front_Stack()){
//			puts("Get!");
			printf("Now stack size %d !\n",StackSize);
			Push_Stack(p=p->lchild);
		}
		puts("Find the left end !");
		printf("Now stack size %d !\n",StackSize);
		Out_Stack();
		puts("out empty p !");
		if(!Empty_Stack()){
			p=Out_Stack();
			puts("View node p !");
//			View_Node(p);
			(*Visit_Node)(p);
			Push_Stack(p=p->rchild);
		}
	}
}
void View_Node(TreeNodePtr n){
	char *ans;
	printf(">>>>>>>>>>>>Now will check node<<<<<<<<<<<<\n :");
	printf("Node :%d Father node :%d ,AttrIndex %d,Gini Index :%f\n",n->NodePosition,n->FNodePos,n->AttrIndex,n->GiniIndex);
	View_Histogram(n->Histogram[1]);
	printf("best split :%f\tacc:%f\tmcc:%f\twacc:%f\twacc:%f\n",n->BestSplit,n->acc,n->mcc,n->wacc,n->wmcc);
	if(n->PureNode)
		printf("Pure Node ?--> NO!\n");
	else
		printf("Pure Node ?--> YES!\n");
}
void BFS_Traverse(Queue *Q,TreeNodePtr root,void ( * Visit_Node)(TreeNodePtr)){
	TreeNodePtr Temp;
	EnQueue(Q,root);
	while(!EmptyQueue(Q)){
		Temp=Get_Front_Queue(Q);
//		View_Node(Temp);
		(*Visit_Node)(Temp);
		if(Temp->lchild) EnQueue(Q,Temp->lchild);
		if(Temp->rchild) EnQueue(Q,Temp->rchild);
	}
}
/*---------------------------------------------------------------------*/

/* Part VI : RI table */

//get the relevance coefficient table
void Initial_RI_Table(){
	int i=0;
	printf("Initial_RI_Table LayerCount %d\n!",LayerCount);
	RI_Table=(RI *)malloc(LayerCount*sizeof(RI));
	for(i=0;i<LayerCount;i++){
		Initial_HIS(RI_Table[i].HIS);
		RI_Table[i].AttrIndex=i;
		RI_Table[i].MCC=0;
		RI_Table[i].MIC=0;
		RI_Table[i].wMIC=0;
		RI_Table[i].ACC=0;
		RI_Table[i].wACC=0;
		RI_Table[i].wMCC=0;
	}
}
void Initial_HIS(int h[2][2]){
	int *p=&h[0][0],i=0;
	while(i<4){
		*(p+i++)=0;
	}
}
void Get_RI_Table(Queue * Q,TreeNodePtr root){
//	Deep_Traverse(root,Cal_Node_AM);
//	Deep_Traverse(root,Cal_Node_AM);
	BFS_Traverse(Q,root,Cal_Node_AM);
	puts("");
	View_RI_Table();
	BFS_Traverse(Q,root,Cal_RI_wAM);
	View_RI_Table();
	Cal_RI_tAM();
	View_RI_Table();
}
void Cal_Node_AM(TreeNodePtr node){
	int layer=node->AttrIndex;
	if(node->PureNode){
		node->mcc=MCC_Cal(node->Histogram[1]);
		node->acc=ACC_Cal(node->Histogram[1]);
		Add_His(RI_Table[layer].HIS,node->Histogram[1]);
	}
}
void Cal_RI_wAM(TreeNodePtr node){
	int layer=node->AttrIndex;
	if(node->PureNode){
		if(node->NodePosition==1){
			node->wmcc=node->mcc;
			node->wacc=node->acc;
		}
		else{
			node->wmcc=(float)Total_HIS(node->Histogram[1])/(float)Total_HIS(RI_Table[layer].HIS)*(node->mcc);
			node->wacc=(float)Total_HIS(node->Histogram[1])/(float)Total_HIS(RI_Table[layer].HIS)*(node->acc);
		}
		RI_Table[layer].wACC+=node->wacc;
		RI_Table[layer].wMCC+=node->wmcc;
	}
}
void Cal_RI_tAM(){
	int i=0;
	for(i=0;i<LayerCount;i++){
		RI_Table[i].ACC=ACC_Cal(RI_Table[i].HIS);
		RI_Table[i].MCC=MCC_Cal(RI_Table[i].HIS);
	}
}
int Total_HIS(int h[2][2]){
	int *p=&h[0][0],i=0,T=0;
	while(i<4){
		T+=*(p+i++);
	}
	return T;
}
void Add_His(int a[2][2],int b[2][2]){	//add b to a and then a is the total of the origianl
	int *p=&a[0][0],*q=&b[0][0],i=0;
	while(i<4){
		*(p+i)=*(p+i)+(*(q+i));
		++i;
	}
}
float MCC_Cal(int a[2][2]){	//only calculate a single histogram
	int TP=a[0][0],FP=a[0][1],FN=a[1][0],TN=a[1][1];
	return (float)(TP*TN-FN*FP)/sqrt((float)(TP+FN)*(TP+FP)*(TN+FN)*(TN+FP));	
}
float ACC_Cal(int a[2][2]){	//only compute a node's condition
	int TP=a[0][0],FP=a[0][1],FN=a[1][0],TN=a[1][1];
	return (float)(TP+TN)/(float)(TP+TN+FP+FN);
}
void View_RI_Table(){
	int i=0;
	printf("\n>>>>>>>>>>>>RI Table<<<<<<<<<<<<<<\n");
	for(i=0;i<LayerCount;i++){
		printf(" %d >>>> MCC:%f\tACC:%f\twMCC:%f\twACC:%f\n",i,RI_Table[i].MCC,\
		RI_Table[i].ACC,RI_Table[i].wMCC,RI_Table[i].wACC);
		View_Histogram(RI_Table[i].HIS);
	}
}
void Destroy_RI_Table(){
	free(RI_Table);
}
/*---------------------------------------------------------------------*/
