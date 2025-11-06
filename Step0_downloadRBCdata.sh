# 定义数据存储路径
cd /Users/lizheng/Desktop/RBC/PNC/PNC_FreeSurfer/freesurfer

DATA_PATH="/Users/lizheng/Desktop/RBC/PNC/PNC_FreeSurfer/freesurfer"  # 替换为你的实际路径
OUTPUT_PATH="/Users/lizheng/Desktop"          # 替换为输出文件的路径
MISSING_FILE_REPORT="${OUTPUT_PATH}/missing_files_report.xlsx"

# 创建输出目录
mkdir -p "$OUTPUT_PATH"

# 初始化缺失文件列表
MISSING_FILES=()

# 定义文件名模式变量
FILE_SUFFIX="brainmeasures"

# 遍历每个子文件夹
for SUBJECT_DIR in "$DATA_PATH"/sub-*; do
    if [ -d "$SUBJECT_DIR" ]; then
        SUBJECT_ID=$(basename "$SUBJECT_DIR")
        TARGET_FILE=$(find "$SUBJECT_DIR" -maxdepth 1 -name "${SUBJECT_ID}_${FILE_SUFFIX}.tsv")

        if [ -z "$TARGET_FILE" ]; then
            # 如果文件不存在，记录缺失的文件夹
            MISSING_FILES+=("$SUBJECT_ID")
        else
            # 下载目标文件
            datalad get "$TARGET_FILE"
        fi
    fi
done

# 如果有缺失文件，生成报告
if [ ${#MISSING_FILES[@]} -gt 0 ]; then
    echo "以下文件夹缺少 *_regionsurfacestats.tsv 文件：" > "${OUTPUT_PATH}/missing_files_report.txt"
    printf "%s\n" "${MISSING_FILES[@]}" >> "${OUTPUT_PATH}/missing_files_report.txt"
    # 使用 Python 将报告转换为 Excel
    python3 <<EOF
import pandas as pd
missing_files = ${MISSING_FILES[@]}
df = pd.DataFrame(missing_files, columns=["Missing Subject ID"])
df.to_excel("${MISSING_FILE_REPORT}", index=False)
EOF
    echo "缺失文件报告已保存到 ${MISSING_FILE_REPORT}"
else
    echo "无缺失文件，所有数据已完整下载。"
fi
