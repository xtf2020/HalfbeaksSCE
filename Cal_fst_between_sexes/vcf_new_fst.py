#!/usr/bin/env python3
"""
直接从VCF文件计算Male vs Female FST (窗口平均值)
用法: python3 vcf_to_fst.py input.vcf pop_file.txt [--per-site]
"""

import sys
import os
import gzip
import logging
import argparse
from collections import defaultdict
import numpy as np

def setup_logging():
    """设置日志记录"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('fst_calculation.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def read_pop_file(pop_file):
    """读取群体文件，返回个体到群体的映射"""
    pop_dict = {}
    male_ids = []
    female_ids = []
    
    with open(pop_file, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 2:
                    indv, pop = parts[0], parts[1].lower()
                    pop_dict[indv] = pop
                    if pop == 'male':
                        male_ids.append(indv)
                    elif pop == 'female':
                        female_ids.append(indv)
    
    logging.info(f"读取群体信息: {len(male_ids)} 雄性, {len(female_ids)} 雌性")
    return pop_dict, male_ids, female_ids

def open_vcf_file(vcf_file):
    """智能打开VCF文件，支持.gz和普通格式"""
    if vcf_file.endswith('.gz'):
        return gzip.open(vcf_file, 'rt')
    else:
        return open(vcf_file, 'r')

def parse_vcf(vcf_file, male_ids, female_ids, maf_threshold=0.2, max_missing=0.9):
    """解析VCF文件，计算等位基因频率
    max_missing: 允许的最大缺失比例，默认0.9表示至少90%个体有数据
    """
    logging.info(f"开始解析VCF文件: {vcf_file}")
    
    # 用于存储每个SNP的信息
    snp_data = []
    
    with open_vcf_file(vcf_file) as f:
        for line in f:
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    # 解析头部，获取个体顺序
                    header = line.strip().split('\t')
                    samples = header[9:]
                    male_indices = [i for i, s in enumerate(samples) if s in male_ids]
                    female_indices = [i for i, s in enumerate(samples) if s in female_ids]
                    total_individuals = len(male_indices) + len(female_indices)
                    required_individuals = int(total_individuals * max_missing)
                    logging.info(f"找到样本索引: 雄性{len(male_indices)}, 雌性{len(female_indices)}")
                    logging.info(f"缺失过滤: 需要至少 {required_individuals}/{total_individuals} 个体有数据")
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            chrom, pos, id_, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            
            # 只处理双等位基因SNP
            if ',' in alt:
                continue
            
            # 解析基因型
            formats = parts[8].split(':')
            gt_index = formats.index('GT') if 'GT' in formats else 0
            
            male_alleles = []
            female_alleles = []
            missing_count = 0  # 新增：统计缺失个体数
            total_counted_individuals = 0  # 新增：统计有数据的个体数
            
            for i, genotype in enumerate(parts[9:]):
                gt_parts = genotype.split(':')
                if len(gt_parts) <= gt_index:
                    missing_count += 1
                    continue
                    
                gt = gt_parts[gt_index]
                # 检查缺失数据
                if '.' in gt or './.' in gt or '.|.' in gt:
                    missing_count += 1
                    continue
                
                # 解析基因型，如 "0/0", "0/1", "1/1"
                alleles = []
                valid_genotype = True
                for a in gt.replace('|', '/').split('/'):
                    if a == '0':
                        alleles.append(0)  # 用0表示ref等位基因
                    elif a == '1':
                        alleles.append(1)  # 用1表示alt等位基因
                    else:
                        valid_genotype = False
                        break
                
                if not valid_genotype:
                    missing_count += 1
                    continue
                
                # 成功解析基因型
                total_counted_individuals += 1
                if i in male_indices:
                    male_alleles.extend(alleles)
                elif i in female_indices:
                    female_alleles.extend(alleles)
            
            # 新增：检查缺失比例
            total_expected_individuals = len(male_indices) + len(female_indices)
            missing_proportion = missing_count / total_expected_individuals if total_expected_individuals > 0 else 1.0
            
            if missing_proportion > (1 - max_missing):
                continue
                
            # 检查样本数量是否正确（基于有数据的个体）
            expected_male_alleles = len([i for i in male_indices if i < len(parts[9:]) and 
                                       not ('.' in parts[9+i].split(':')[gt_index] if gt_index < len(parts[9+i].split(':')) else True)]) * 2
            expected_female_alleles = len([i for i in female_indices if i < len(parts[9:]) and 
                                         not ('.' in parts[9+i].split(':')[gt_index] if gt_index < len(parts[9+i].split(':')) else True)]) * 2
            
            # 计算MAF并过滤
            all_alleles = male_alleles + female_alleles
            if len(all_alleles) == 0:
                continue
                
            alt_freq = sum(all_alleles) / len(all_alleles)
            maf = min(alt_freq, 1 - alt_freq)
            
            if maf < maf_threshold:
                continue
            
            snp_data.append({
                'chrom': chrom,
                'pos': int(pos),
                'male_alleles': male_alleles,
                'female_alleles': female_alleles,
                'maf': maf,
                'missing_proportion': missing_proportion  # 新增：记录缺失比例
            })
    
    logging.info(f"成功解析 {len(snp_data)} 个SNP (MAF >= {maf_threshold}, 缺失比例 <= {1-max_missing:.1%})")
    return snp_data

def calculate_weir_cockerham_fst(male_alleles, female_alleles):
    """计算真实的Weir & Cockerham FST估计值 - 修正版"""
    n_m = len(male_alleles)  # 雄性等位基因数
    n_f = len(female_alleles)  # 雌性等位基因数
    
    if n_m == 0 or n_f == 0:
        return np.nan
    
    # 等位基因频率
    p_m = np.mean(male_alleles)
    p_f = np.mean(female_alleles)
    
    # Weir & Cockerham 计算 - 关键修正！
    r = 2  # 群体数
    N = n_m + n_f  # 总等位基因数
    
    # 加权平均频率
    p_bar = (n_m * p_m + n_f * p_f) / N
    
    # 平方和
    SSG = n_m * (p_m - p_bar)**2 + n_f * (p_f - p_bar)**2
    SSI = n_m * p_m * (1 - p_m) + n_f * p_f * (1 - p_f)
    
    # 均方 - 关键修正！
    MSG = SSG / (r - 1)
    MSI = SSI / (N - r)  # 分母是 N - r，不是 sum(n_i-1)
    
    # 校正样本量
    n_c = (N - (n_m**2 + n_f**2) / N) / (r - 1)
    
    if n_c <= 0 or MSI <= 0:
        return np.nan
    
    # 方差分量
    sigma_a2 = (MSG - MSI) / n_c
    sigma_b2 = MSI
    
    # 处理负值
    if sigma_a2 < 0:
        sigma_a2 = 0
    
    if sigma_a2 + sigma_b2 == 0:
        return 0.0
    
    # Weir & Cockerham FST
    theta = sigma_a2 / (sigma_a2 + sigma_b2)
    return max(0.0, theta)  # 负值截断为0

def calculate_nei_gst(male_alleles, female_alleles):
    """计算Nei's G_ST"""
    p_m = np.mean(male_alleles)
    p_f = np.mean(female_alleles)
    
    # 预期杂合度
    H_S_m = 2 * p_m * (1 - p_m)
    H_S_f = 2 * p_f * (1 - p_f)
    H_S = (H_S_m + H_S_f) / 2
    
    # 总群体杂合度
    # all_alleles = male_alleles + female_alleles
    # p_total = np.mean(all_alleles)
    # H_T = 2 * p_total * (1 - p_total)

    # 总群体杂合度
    H_T_x = p_m + p_f
    H_T_y = (1 - p_m) + (1 - p_f)
    H_T = (H_T_x * H_T_y) / 2
    
    if H_T == 0:
        return 0.0
    
    gst = (H_T - H_S) / H_T
    return max(0.0, gst)  # 负值截断为0

def calculate_hudson_fst(male_alleles, female_alleles):
    """计算Hudson FST"""
    p_m = np.mean(male_alleles)
    p_f = np.mean(female_alleles)
    
    H_W = 0.5 * (2 * p_f * (1 - p_f) + 2 * p_m * (1 - p_m))
    H_B = p_f * (1 - p_m) + p_m * (1 - p_f)
    
    if H_B == 0:
        return 0.0
    
    fst = 1 - (H_W / H_B)
    return max(0.0, fst)  # 负值截断为0

def calculate_fst_window(snps_in_window, method='weir'):
    """计算窗口内SNP的平均FST"""
    fst_values = []
    
    for snp in snps_in_window:
        male_alleles = snp['male_alleles']
        female_alleles = snp['female_alleles']
        
        if method == 'weir':
            fst = calculate_weir_cockerham_fst(male_alleles, female_alleles)
        elif method == 'nei':
            fst = calculate_nei_gst(male_alleles, female_alleles)
        elif method == 'hudson':
            fst = calculate_hudson_fst(male_alleles, female_alleles)
        else:
            fst = np.nan
        
        if not np.isnan(fst):
            fst_values.append(fst)
    
    if not fst_values:
        return np.nan
    
    return np.mean(fst_values)

def calculate_expected_fst(n_male, n_female):
    """计算在当前样本量下的预期FST值 - 修正版"""
    logging.info("计算预期FST值...")
    
    # 模拟完全固定的SNP：雌性A/A，雄性A/C
    n_f = n_female * 2
    n_m = n_male * 2
    
    # 创建模拟数据 - 修正：雄性为交替的0和1，模拟A/C杂合
    male_alleles_fixed = [0, 1] * n_male  # A/C 杂合
    female_alleles_fixed = [0] * n_f      # A/A 纯合
    
    expected_weir = calculate_weir_cockerham_fst(male_alleles_fixed, female_alleles_fixed)
    expected_nei = calculate_nei_gst(male_alleles_fixed, female_alleles_fixed)
    expected_hudson = calculate_hudson_fst(male_alleles_fixed, female_alleles_fixed)
    
    return {
        'weir': expected_weir,
        'nei': expected_nei, 
        'hudson': expected_hudson
    }

def main():
    parser = argparse.ArgumentParser(description='从VCF计算Male vs Female FST')
    parser.add_argument('vcf_file', help='输入VCF文件(.vcf或.vcf.gz)')
    parser.add_argument('pop_file', help='群体文件')
    parser.add_argument('--window-size', type=int, default=50000, help='窗口大小(bp)，默认50kb')
    parser.add_argument('--per-site', action='store_true', help='输出每个SNP的FST值')
    parser.add_argument('--maf', type=float, default=0.2, help='MAF阈值，默认0.2')
    parser.add_argument('--max-missing', type=float, default=0.9, 
                       help='最大缺失比例，默认0.9表示至少90%%个体有数据，1.0表示不进行缺失过滤')
    
    args = parser.parse_args()
    
    # 设置日志
    setup_logging()
    
    # 读取群体文件
    pop_dict, male_ids, female_ids = read_pop_file(args.pop_file)
    
    if not male_ids or not female_ids:
        logging.error("群体文件中必须包含male和female个体")
        sys.exit(1)
    
    # 计算预期FST值
    expected = calculate_expected_fst(len(male_ids), len(female_ids))
    
    logging.info("=== 在当前样本量下的预期FST值 (完全固定SNP) ===")
    logging.info(f"样本量: {len(male_ids)} 雄性, {len(female_ids)} 雌性")
    logging.info(f"过滤参数: MAF >= {args.maf}, 缺失比例 <= {1-args.max_missing:.1%}")
    logging.info(f"Weir & Cockerham FST: {expected['weir']:.6f}")
    logging.info(f"Nei's G_ST: {expected['nei']:.6f}")
    logging.info(f"Hudson FST: {expected['hudson']:.6f}")
    
    # 验证已知情况
    test_cases = [(20, 20), (21, 27)]
    logging.info("验证计算 (应与vcftools匹配):")
    for n_m, n_f in test_cases:
        test_fst = calculate_expected_fst(n_m, n_f)
        logging.info(f"  {n_m}雄, {n_f}雌: Weir FST = {test_fst['weir']:.6f}")
    
    logging.info("=" * 50)
    
    # 解析VCF文件 - 传入max_missing参数
    snp_data = parse_vcf(args.vcf_file, male_ids, female_ids, args.maf, args.max_missing)
    
    if not snp_data:
        logging.error("没有找到可用的SNP数据")
        sys.exit(1)
    
    # 按染色体和位置排序
    snp_data.sort(key=lambda x: (x['chrom'], x['pos']))
    
    # 计算FST
    if args.per_site:
        # 输出每个SNP的FST
        with open('fst_per_site.txt', 'w') as f:
            f.write("CHROM\tPOS\tFST_WEIR\tFST_NEI\tFST_HUDSON\tMAF\tMISSING_PROP\n")
            for snp in snp_data:
                fst_weir = calculate_weir_cockerham_fst(snp['male_alleles'], snp['female_alleles'])
                fst_nei = calculate_nei_gst(snp['male_alleles'], snp['female_alleles'])
                fst_hudson = calculate_hudson_fst(snp['male_alleles'], snp['female_alleles'])
                
                f.write(f"{snp['chrom']}\t{snp['pos']}\t{fst_weir:.6f}\t{fst_nei:.6f}\t{fst_hudson:.6f}\t{snp['maf']:.4f}\t{snp.get('missing_proportion', 0):.4f}\n")
        
        logging.info("每个SNP的FST结果已保存到: fst_per_site.txt")
    
    else:
        # 计算窗口FST
        window_size = args.window_size
        logging.info(f"计算 {window_size/1000}kb 窗口FST")
        
        with open('fst_window_results.txt', 'w') as f:
            f.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tFST_WEIR\tFST_NEI\tFST_HUDSON\n")
            
            # 按染色体分组
            chromosomes = set(snp['chrom'] for snp in snp_data)
            
            for chrom in sorted(chromosomes):
                chrom_snps = [s for s in snp_data if s['chrom'] == chrom]
                if not chrom_snps:
                    continue
                
                max_pos = max(s['pos'] for s in chrom_snps)
                
                # vcftools风格的滑动窗口：从1开始，连续窗口
                first_window_start = 1
                max_pos_rounded = ((max_pos + window_size - 1) // window_size) * window_size
                
                for start in range(first_window_start, max_pos_rounded + 1, window_size):
                    end = start + window_size - 1
                    snps_in_window = [s for s in chrom_snps if start <= s['pos'] <= end]
                    
                    if len(snps_in_window) < 3:  # 至少3个SNP
                        continue
                    
                    fst_weir = calculate_fst_window(snps_in_window, 'weir')
                    fst_nei = calculate_fst_window(snps_in_window, 'nei')
                    fst_hudson = calculate_fst_window(snps_in_window, 'hudson')
                    
                    if not np.isnan(fst_weir):
                        f.write(f"{chrom}\t{start}\t{end}\t{len(snps_in_window)}\t{fst_weir:.6f}\t{fst_nei:.6f}\t{fst_hudson:.6f}\n")
        
        logging.info(f"窗口FST结果已保存到: fst_window_results.txt")
    
    logging.info("分析完成!")

if __name__ == "__main__":
    main()