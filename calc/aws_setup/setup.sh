#!/usr/bin/env bash
# =============================================================================
# ONE-TIME AWS setup — run with YOUR ADMIN creds on YOUR machine (not in the
# sandbox). Creates a least-privilege IAM user + EC2 instance role + private
# encrypted results bucket for the alpha-O job, and prints the scoped access key
# to paste into the Claude Code environment config.
#
# Policies are adversarially verified: sufficient (no AccessDenied) AND minimal.
# No SSH / no key-pair / no inbound; plain AES256 bucket (no KMS perms needed).
#
#   REGION=us-east-1 ./setup.sh
# Override any of: REGION BUCKET PREFIX ROLE_NAME PROFILE_NAME USER_NAME TAG
# =============================================================================
set -euo pipefail

REGION="${REGION:-us-east-1}"
ACCOUNT_ID="$(aws sts get-caller-identity --query Account --output text)"
BUCKET="${BUCKET:-alpha-o-results-$ACCOUNT_ID}"
PREFIX="${PREFIX:-alpha-o}"
ROLE_NAME="${ROLE_NAME:-alpha-o-instance-role}"
PROFILE_NAME="${PROFILE_NAME:-alpha-o-instance-profile}"
USER_NAME="${USER_NAME:-alpha-o-agent}"
TAG="${TAG:-sandbox-agent}"
echo "region=$REGION account=$ACCOUNT_ID bucket=$BUCKET prefix=$PREFIX role=$ROLE_NAME user=$USER_NAME"

render() { sed -e "s/REGION/$REGION/g" -e "s/ACCOUNT_ID/$ACCOUNT_ID/g" \
               -e "s/ROLE_NAME/$ROLE_NAME/g" -e "s/BUCKET/$BUCKET/g" \
               -e "s/PREFIX/$PREFIX/g" -e "s/sandbox-agent/$TAG/g" "$1" | grep -v '"_comment"'; }
HERE="$(cd "$(dirname "$0")" && pwd)"

# --- private, encrypted results bucket ---
if ! aws s3api head-bucket --bucket "$BUCKET" 2>/dev/null; then
  if [ "$REGION" = us-east-1 ]; then
    aws s3api create-bucket --bucket "$BUCKET" --region "$REGION"
  else
    aws s3api create-bucket --bucket "$BUCKET" --region "$REGION" \
      --create-bucket-configuration LocationConstraint="$REGION"
  fi
fi
aws s3api put-public-access-block --bucket "$BUCKET" --public-access-block-configuration \
  BlockPublicAcls=true,IgnorePublicAcls=true,BlockPublicPolicy=true,RestrictPublicBuckets=true
aws s3api put-bucket-encryption --bucket "$BUCKET" --server-side-encryption-configuration \
  '{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'

# --- EC2 instance role + profile (S3 write to one prefix only) ---
aws iam create-role --role-name "$ROLE_NAME" --assume-role-policy-document \
  '{"Version":"2012-10-17","Statement":[{"Effect":"Allow","Principal":{"Service":"ec2.amazonaws.com"},"Action":"sts:AssumeRole"}]}' 2>/dev/null || true
aws iam put-role-policy --role-name "$ROLE_NAME" --policy-name s3-write \
  --policy-document "$(render "$HERE/instance_role_policy.json")"
aws iam create-instance-profile --instance-profile-name "$PROFILE_NAME" 2>/dev/null || true
aws iam add-role-to-instance-profile --instance-profile-name "$PROFILE_NAME" --role-name "$ROLE_NAME" 2>/dev/null || true

# --- scoped IAM user + inline policy ---
aws iam create-user --user-name "$USER_NAME" 2>/dev/null || true
aws iam put-user-policy --user-name "$USER_NAME" --policy-name alpha-o-run \
  --policy-document "$(render "$HERE/user_policy.json")"

echo; echo "=================================================================="
echo "ACCESS KEY — paste into the Claude Code environment config, then DELETE"
echo "this key after the run (aws iam delete-access-key --user-name $USER_NAME --access-key-id <id>)"
echo "=================================================================="
read AKID SECRET < <(aws iam create-access-key --user-name "$USER_NAME" \
  --query 'AccessKey.[AccessKeyId,SecretAccessKey]' --output text)
cat <<ENV

# --- Environment variables (paste into the env config, .env format, no quotes) ---
AWS_ACCESS_KEY_ID=$AKID
AWS_SECRET_ACCESS_KEY=$SECRET
AWS_DEFAULT_REGION=$REGION
ALPHA_O_S3=s3://$BUCKET/$PREFIX
ALPHA_O_INSTANCE_PROFILE=$PROFILE_NAME
ALPHA_O_TAG=$TAG
ENV

echo
echo "RECOMMENDED: a budget guardrail, e.g.:"
echo "  aws budgets create-budget --account-id $ACCOUNT_ID --budget '{\"BudgetName\":\"alpha-o\",\"BudgetLimit\":{\"Amount\":\"50\",\"Unit\":\"USD\"},\"TimeUnit\":\"MONTHLY\",\"BudgetType\":\"COST\"}'"
